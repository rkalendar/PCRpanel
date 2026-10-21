// Identifies repeated k-mers across a set of DNA sequences. 
public final class MaskingSequences {

    /**
     * Minimum occurrence count for a k-mer to be considered a repeat.
     */
    public static final int MIN_COUNT = 2;

    /**
     * Sentinel index used by {@code tables.dx2} for ambiguous bases (e.g. 'N').
     */
    private static final int AMBIGUOUS_BASE = 4;
    
 
    /**
     * @param sequences array of DNA sequence strings
     * @param kmerSize length of k-mers to count
     * @return LongIntHashMap где value >= MIN_COUNT  
     */
    public LongIntHashMap mask(String[] sequences, int kmerSize) {
        long totalBases = 0L;
        for (String seq : sequences) {
            totalBases += seq.length();
        }
        // Initial capacity — the nearest power of two, no more than 4M
        int initialCapacity = (int) Math.min(totalBases, 1 << 22);
        LongIntHashMap counts = new LongIntHashMap(initialCapacity);

        for (String seq : sequences) {
            countKmers(counts, seq, kmerSize);
        }
        return counts;
    }

    /**
     * Counts all valid k-mers (no ambiguous bases) on both strands of a single
     * sequence, merging counts into {@code counts}.
     */
    private void countKmers(LongIntHashMap counts, String seq, int kmerSize) {
        byte[] forward = toEncodedBases(seq);
        int length = forward.length;

        countKmersOnStrand(counts, forward, length, kmerSize);

        forward = toEncodedBases(dna.ComplementDNA(seq));
        countKmersOnStrand(counts, forward, length, kmerSize);
    }

    /**
     * Sliding-window k-mer counting on a single encoded strand.
     */
    private void countKmersOnStrand(LongIntHashMap counts, byte[] bases,
            int length, int kmerSize) {
        int ambiguousCount = 0;

        for (int i = 0; i < kmerSize - 1; i++) {
            if (bases[i] == AMBIGUOUS_BASE) {
                ambiguousCount++;
            }
        }

        for (int i = kmerSize - 1; i < length; i++) {
            if (bases[i] == AMBIGUOUS_BASE) {
                ambiguousCount++;
            }

            if (ambiguousCount == 0) {
                long code = encodeKmer(bases, i - kmerSize + 1, kmerSize);
                int current = counts.get(code);
                if (current == LongIntHashMap.NO_VALUE) {
                    counts.put(code, 1);
                } else if (current < MIN_COUNT) {
                    counts.put(code, current + 1);
                }
            }

            if (bases[i + 1 - kmerSize] == AMBIGUOUS_BASE) {
                ambiguousCount--;
            }
        }
    }

    /**
     * Converts a DNA string to an encoded byte array using the lookup table
     * {@code tables.dx2} (expected mapping: A→0, C→1, G→2, T→3, other→4).
     */
    private static byte[] toEncodedBases(String seq) {
        byte[] raw = seq.getBytes();
        for (int i = 0; i < raw.length; i++) {
            raw[i] = tables.dx2[raw[i]];
        }
        return raw;
    }

    /**
     * Base-5 k-mer encoding. Collision-free for the DNA+gap alphabet; stays
     * positive for k <= 27.
     */
    private long encodeKmer(byte[] bases, int start, int kmerSize) {
        long code = 0L;
        for (int i = start; i < start + kmerSize; i++) {
            code = code * 5 + bases[i];
        }
        return code;
    }
}
