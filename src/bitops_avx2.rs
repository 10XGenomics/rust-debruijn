use std::arch::x86_64::*;

/// pack the lowest 2 bits of each byte of a 32 byte slice into a u64
/// all bytes must be have values in the range 0..3 incorrect results will be returned.
/// The first byte in the m256 will become the highest 2bits in the output, consistent
/// with the lexicographic sorting convention in this crate.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn pack_32_bases(bases: __m256i) -> u64 {
    // bases = d c b a

    let reverse_mask = _mm256_set_epi8(
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15,
    );

    // step 1: reverse all bytes within lanes
    // reversed = c d a b
    let reversed = _mm256_shuffle_epi8(bases, reverse_mask);

    // step 2: use a lane crossing permute to reverse lanes and
    // swap the middle two 64-bit chunks
    // permuted = a c b d
    let permuted = _mm256_permute4x64_epi64(reversed, 0b01_11_00_10);

    // step 3: interleave the bytes that contain the first and second bits
    let first_bits = _mm256_slli_epi16(permuted, 7);
    let second_bits = _mm256_slli_epi16(permuted, 6);

    // i(x) = interleave first and second bits of each byte of x
    // lo_half = i(c) i(d)
    let lo_half = _mm256_unpacklo_epi8(first_bits, second_bits);

    // hi_half = i(a) i(b)
    let hi_half = _mm256_unpackhi_epi8(first_bits, second_bits);

    // step 4: extract bits using movemask (zero extend)
    let packed_lo = (_mm256_movemask_epi8(lo_half) as u32) as u64;
    let packed_hi = (_mm256_movemask_epi8(hi_half) as u32) as u64;

    (packed_hi << 32) | packed_lo
}

/// Convert a slice of 32 ACGT bytes into a __m256i using 0-4 encoding.
/// Second element in tuple will be false if any of the bases are not in [aAcCgGtT].
/// Bases not in [aAcCgGtT] will be converted to A / 0.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn convert_bases(bytes: &[u8]) -> (__m256i, bool) {
    assert!(bytes.len() == 32);

    // a lookup table to map a number from 0..128 to a single bit
    let hi_lut = {
        let mut lut_hi = 0i64;
        lut_hi |= 1i64 << ((b'A' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'C' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'G' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'T' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'a' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'c' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'g' as i64) - 64i64);
        lut_hi |= 1i64 << ((b't' as i64) - 64i64);
        _mm256_set_epi64x(lut_hi, 0i64, lut_hi, 0i64)
    };

    let lo_lut = _mm256_set_epi8(
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
    );

    let lo_mask = _mm256_set1_epi8(0b00001111);

    // convert ACGT (case-insensitive) to a 2-bit representation by looking up low 4 bits
    let lut = _mm256_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0, /* G */ 2, 0, 0, /* T */ 3, /* C */ 1, 0,
        /* A */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* G */ 2, 0, 0, /* T */ 3,
        /* C */ 1, 0, /* A */ 0, 0,
    );

    let input = _mm256_loadu_si256(bytes.as_ptr() as *const __m256i);

    // step 1: lookup a byte using the high 4 bits of each character
    let hi = _mm256_and_si256(_mm256_srli_epi16(input, 3), lo_mask);
    let hi_lookup = _mm256_shuffle_epi8(hi_lut, hi);

    // step 2: convert the low 3 bits of each character to a mask
    // byte x is converted to (1 << x)
    let lo_lookup = _mm256_shuffle_epi8(lo_lut, input);

    // step 3: AND the byte derived from the high 4 bits and the mask derived from the low 3 bits
    let mask = _mm256_cmpeq_epi8(
        _mm256_and_si256(lo_lookup, hi_lookup),
        _mm256_setzero_si256(),
    );
    let valid = _mm256_testc_si256(_mm256_setzero_si256(), mask) != 0;

    // step 4: use lookup table to convert nucleotides to their 2-bit representation,
    // while zeroing out invalid characters (shuffle returns a zero byte if MSB of a byte is 1)
    let shuffled = _mm256_shuffle_epi8(lut, input);
    let res = _mm256_andnot_si256(mask, shuffled);

    (res, valid)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{base_to_bits, bits_to_ascii};
    use crate::{test, Kmer};

    #[test]
    fn test_convert_bytes_for_debug() {
        let cc = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let gg = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
        let tt = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let at = b"ATATATATATATATATATATATATATATATAT";
        let aatt = b"AAAAAAAAAAAAAAAATTTTTTTTTTTTTTTT";
        let acgt = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let acgt8 = b"AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT";
        let acgt4 = b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT";

        for dna_ascii in vec![cc, gg, tt, at, aatt, acgt, acgt4, acgt8] {
            //let dna_bytes = test::random_dna(32);
            //let dna_ascii: Vec<_> = dna_bytes.iter().map(|bits| bits_to_ascii(*bits)).collect();
            let dna_bytes: Vec<_> = dna_ascii.iter().map(|base| base_to_bits(*base)).collect();

            let dna_str = std::str::from_utf8(dna_ascii).unwrap();
            let (simd_bytes, valid) = unsafe { convert_bases(dna_ascii) };

            let packed = unsafe { pack_32_bases(simd_bytes) };
            let kmer = crate::kmer::Kmer32::from_u64(packed);
            println!("orig: {} \nresult: {:?}", dna_str, kmer);

            assert!(valid);

            for i in 0..32 {
                let v = extract_byte(simd_bytes, i);
                assert_eq!(v, dna_bytes[i]);
            }
        }
    }

    #[test]
    fn test_convert_bytes_random() {
        for _ in 0..1000 {
            let dna_bytes = test::random_dna(32);
            let dna_ascii: Vec<_> = dna_bytes.iter().map(|bits| bits_to_ascii(*bits)).collect();

            let dna_str = std::str::from_utf8(&dna_ascii).unwrap();
            let (simd_bytes, valid) = unsafe { convert_bases(&dna_ascii) };

            let packed = unsafe { pack_32_bases(simd_bytes) };
            let kmer = crate::kmer::Kmer32::from_u64(packed);
            println!("orig: {} \nresult: {:?}", dna_str, kmer);

            assert_eq!(dna_str, &kmer.to_string());
            assert!(valid);

            for i in 0..32 {
                let v = extract_byte(simd_bytes, i);
                assert_eq!(v, dna_bytes[i]);
            }
        }
    }

    #[test]
    fn test_invalid_bases() {
        for i in 0..1000 {
            let dna_bytes = test::random_dna(32);
            let mut dna_ascii: Vec<_> = dna_bytes.iter().map(|bits| bits_to_ascii(*bits)).collect();

            for j in 0..50 {
                let pos = i * j % 32;
                dna_ascii[pos] = (i * j % 256) as u8;

                let (_, simd_valid) = unsafe { convert_bases(&dna_ascii) };

                let true_valid = !dna_ascii
                    .iter()
                    .any(|b| crate::dna_only_base_to_bits(*b).is_none());

                if simd_valid != true_valid {
                    println!("{:?}", dna_ascii);
                }
                assert_eq!(simd_valid, true_valid);
            }
        }
    }

    fn _get(v: __m256i) -> i64 {
        unsafe { _mm256_extract_epi64(v, 0) }
    }

    fn extract_byte(v: __m256i, idx: usize) -> u8 {
        let idx = idx as i8;

        let res = unsafe {
            match idx {
                0 => _mm256_extract_epi8(v, 0),
                1 => _mm256_extract_epi8(v, 1),
                2 => _mm256_extract_epi8(v, 2),
                3 => _mm256_extract_epi8(v, 3),
                4 => _mm256_extract_epi8(v, 4),
                5 => _mm256_extract_epi8(v, 5),
                6 => _mm256_extract_epi8(v, 6),
                7 => _mm256_extract_epi8(v, 7),
                8 => _mm256_extract_epi8(v, 8),
                9 => _mm256_extract_epi8(v, 9),
                10 => _mm256_extract_epi8(v, 10),
                11 => _mm256_extract_epi8(v, 11),
                12 => _mm256_extract_epi8(v, 12),
                13 => _mm256_extract_epi8(v, 13),
                14 => _mm256_extract_epi8(v, 14),
                15 => _mm256_extract_epi8(v, 15),
                16 => _mm256_extract_epi8(v, 16),
                17 => _mm256_extract_epi8(v, 17),
                18 => _mm256_extract_epi8(v, 18),
                19 => _mm256_extract_epi8(v, 19),
                20 => _mm256_extract_epi8(v, 20),
                21 => _mm256_extract_epi8(v, 21),
                22 => _mm256_extract_epi8(v, 22),
                23 => _mm256_extract_epi8(v, 23),
                24 => _mm256_extract_epi8(v, 24),
                25 => _mm256_extract_epi8(v, 25),
                26 => _mm256_extract_epi8(v, 26),
                27 => _mm256_extract_epi8(v, 27),
                28 => _mm256_extract_epi8(v, 28),
                29 => _mm256_extract_epi8(v, 29),
                30 => _mm256_extract_epi8(v, 30),
                31 => _mm256_extract_epi8(v, 31),
                _ => panic!("bad index"),
            }
        };

        res as u8
    }

    fn _print64b(v: __m256i) -> String {
        unsafe {
            format!(
                "{:#b} {:#b} {:#b} {:#b}",
                _mm256_extract_epi64(v, 0),
                _mm256_extract_epi64(v, 1),
                _mm256_extract_epi64(v, 2),
                _mm256_extract_epi64(v, 3)
            )
        }
    }
}
