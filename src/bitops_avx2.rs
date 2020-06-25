//#[cfg(test)]
//#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
//#[target_feature(enable = "avx2")]
//mod avx2 {

use std::arch::x86_64::*;

fn get(v: __m256i) -> i64 {
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

fn print64b(v: __m256i) -> String {
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

/// pack the lowest 2 bits of each byte of a 32 byte slice into a u64
/// all bytes must be have values in the range 0..3 incorrect results will be returned.
/// The first byte in the m256 will become the highest 2bits in the output, consistent
/// with the lexicographic sorting convention in this crate.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn pack_32_bases(bases: __m256i) -> u64 {
    // mulitplier to move the bottom two bits from 4 bytes together
    // using a single 32bit multiplication.
    let pack_const: i32 = (1 << 6) + (1 << (4 + 8)) + (1 << (2 + 16)) + (1 << 24);
    let pack_const_256 = _mm256_set1_epi32(pack_const as i32);

    // extract 4th byte of each 64-bit chunk. this is the byte in the
    // packing multiplication that contains the 2-bit values all together.
    let byte_sel_half = _mm_set_epi8(
        0 + 3,
        0 + 3,
        8 + 3,
        8 + 3,
        16 + 3,
        16 + 3,
        24 + 3,
        24 + 3,
        0 + 3,
        0 + 3,
        8 + 3,
        8 + 3,
        16 + 3,
        16 + 3,
        24 + 3,
        24 + 3,
    );
    let byte_sel = _mm256_set_m128i(byte_sel_half, byte_sel_half);

    // extract bits from the lower 32bit of each 64.
    let pack_low = _mm256_mul_epu32(bases, pack_const_256);

    // move the top 32 down to the bottom in each 64.
    let bases_shuf = _mm256_shuffle_epi32(bases, 0b11110101);

    // extract bits from the upper 32bit of each 64.
    let pack_high = _mm256_mul_epu32(bases_shuf, pack_const_256);

    let compacted_low = _mm256_shuffle_epi8(pack_low, byte_sel);
    let compacted_high = _mm256_shuffle_epi8(pack_high, byte_sel);

    // blend low and high bytes back together
    let blend_mask = _mm256_set_epi8(
        -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0,
        -1, 0, -1, 0, -1, 0,
    );
    let compacted = _mm256_blendv_epi8(compacted_high, compacted_low, blend_mask);

    // bring the 2 32-bit chunks together
    let joined = _mm256_permutevar8x32_epi32(compacted, _mm256_setr_epi32(4, 0, 7, 7, 7, 7, 7, 7));

    // reverse the byte in groups of 2 bits
    let rev1 = _mm256_or_si256(
        _mm256_slli_epi16(_mm256_and_si256(joined, packed_u16(0x3333)), 2),
        _mm256_and_si256(_mm256_srli_epi16(joined, 2), packed_u16(0x3333)),
    );

    let rev2 = _mm256_or_si256(
        _mm256_slli_epi16(_mm256_and_si256(rev1, packed_u16(0x0F0F)), 4),
        _mm256_and_si256(_mm256_srli_epi16(rev1, 4), packed_u16(0x0F0F)),
    );

    let packed64 = _mm256_extract_epi64(rev2, 0);
    packed64 as u64
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn packed_byte(byte: u8) -> __m256i {
    _mm256_set1_epi8(byte as i8)
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn packed_u16(v: u16) -> __m256i {
    _mm256_set1_epi16(v as i16)
}

/// Covert a slice of 32 ACGT bytes into a __m256i using 0-4 encoding.
/// Second element in tuple will be false if any of the bases are not in [aAcCgGtT].
/// Bases not in [aAcCgGtT] will be converted to A / 0.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn convert_bases(bytes: &[u8]) -> (__m256i, bool) {
    assert!(bytes.len() == 32);

    let input = _mm256_loadu_si256(bytes.as_ptr() as *const __m256i);

    let a1 = _mm256_cmpeq_epi8(input, packed_byte(b'a'));
    let a2 = _mm256_cmpeq_epi8(input, packed_byte(b'A'));
    let a = _mm256_or_si256(a1, a2);

    let c1 = _mm256_cmpeq_epi8(input, packed_byte(b'c'));
    let c2 = _mm256_cmpeq_epi8(input, packed_byte(b'C'));
    let c = _mm256_or_si256(c1, c2);

    let g1 = _mm256_cmpeq_epi8(input, packed_byte(b'g'));
    let g2 = _mm256_cmpeq_epi8(input, packed_byte(b'G'));
    let g = _mm256_or_si256(g1, g2);

    let t1 = _mm256_cmpeq_epi8(input, packed_byte(b't'));
    let t2 = _mm256_cmpeq_epi8(input, packed_byte(b'T'));
    let t = _mm256_or_si256(t1, t2);

    let mut res = packed_byte(0);
    let mut valid = packed_byte(0);

    res = _mm256_blendv_epi8(res, packed_byte(0), a);
    valid = _mm256_or_si256(valid, a);

    res = _mm256_blendv_epi8(res, packed_byte(1), c);
    valid = _mm256_or_si256(valid, c);

    res = _mm256_blendv_epi8(res, packed_byte(2), g);
    valid = _mm256_or_si256(valid, g);

    res = _mm256_blendv_epi8(res, packed_byte(3), t);
    valid = _mm256_or_si256(valid, t);

    let vv = _mm256_testc_si256(valid, packed_byte(255));

    return (res, vv == 1);
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
}
