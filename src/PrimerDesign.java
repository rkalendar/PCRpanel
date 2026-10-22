import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Designs PCR primers for given DNA sequences, applying masking against a
 * reference genome and optional homology filtering. Exon regions are merged
 * when they fall within PCR size constraints, and candidate forward/reverse
 * primers are evaluated for melting temperature, linguistic complexity,
 * self-dimerization, and simple sequence repeat content.
 */
public class PrimerDesign {

    // ── Constants ────────────────────────────────────────────────────────
    /** k-mer size for de novo repeat detection. Configurable via 'repeatKmer='. */
    public static int GENOME_KMER_SIZE = 18;
    /** Repeat-masked positions tolerated inside a candidate primer.
     *  Configurable via 'maskedBases='; 0 reproduces the original behaviour. */
    public static int MASK_TOLERANCE = 0;
    private static final int MAX_EXTRA_SEQ = 500;
    private static final int RELAXED_LC_OFFSET = 0;// 5;
    private static final int RELAXED_MAXLEN_OFFSET = 0; // 5;
    private static final int MIN_PRIMERS_BEFORE_RELAXING = 0;//2;

    /**
     * Simple Sequence Repeat patterns that disqualify a primer candidate.
     */
    private static final Set<String> SSR_PATTERNS = Set.of("gggg", "cccc", "aaaa", "tttt", "cgcgc", "gcgcg", "tatat", "atata", "tgtgt", "gtgtg", "cacaca", "acacac", "tctctc", "ctctct", "agagag", "gagaga");

    // ── Configuration ────────────────────────────────────────────────────
    private final boolean useHomologyFilter;
    private final String sequence;
    private final String sequenceName;
    private final String referenceSeq;
    private final String fastaSeq;
    private final int minLinguisticComplexity;
    private final int minPrimerLength;
    private final int maxPrimerLength;
    private final int minMeltingTemp;
    private final int maxMeltingTemp;

    private final int[] exonBoundaries;
    private final String[] existingPrimers;
    private final LongIntHashMap genomeKmerSet;

    // ── Results ──────────────────────────────────────────────────────────
    private PrimersCollector[] forwardPrimers;
    private PrimersCollector[] reversePrimers;
    private int[] widenedFlanks; // per region: search flank (bp) used by the rescue, 0 = default
    private List<Pair> allPairs = new ArrayList<>(); // every valid primer combination, for closing gaps
    private SearchMap search = new SearchMap(false);

    /**
     * What the primer search found along the target sequence, for the panel
     * map: where no primer can bind and why, where the primers that pass the
     * filters lie, and why a region got no primer pair. Coordinates are 0-based
     * half-open on the sense strand.
     */
    public static final class SearchMap {

        /** k-mers that occur more than once in the target, the reference sequences or the other targets */
        public final List<int[]> repeats = new ArrayList<>();
        /** k-mers that occur more than once in the genome reference */
        public final List<int[]> genome = new ArrayList<>();
        /**
         * Bases that no primer of the minimum length can include: its linguistic
         * complexity is below minLC, or it holds a simple sequence repeat.
         * Checked within maxPCR of the targets only.
         */
        public final List<int[]> lowComplexity = new ArrayList<>();
        /** false when no genome reference was given, so {@link #genome} is empty by default */
        public final boolean genomeChecked;
        /** the candidate primers of every design region */
        public final List<Site> sites = new ArrayList<>();
        /** per design region, why it has no primer pair; null where it has one */
        public final List<String> pairing = new ArrayList<>();

        SearchMap(boolean genomeChecked) {
            this.genomeChecked = genomeChecked;
        }

        /** Why design region {@code region} has no primer pair, or null. */
        public String pairing(int region) {
            return (region >= 0 && region < pairing.size()) ? pairing.get(region) : null;
        }
    }

    /** A candidate primer: one that passed every filter, paired or not. */
    public static final class Site {

        public final boolean forward;
        public final int start;     // binding site on the sense strand, 0-based half-open
        public final int end;
        public final String name;
        public final double tm;
        public final int region;    // design region, 0-based as Pair.exonid

        Site(boolean forward, int start, int end, String name, double tm, int region) {
            this.forward = forward;
            this.start = start;
            this.end = end;
            this.name = name;
            this.tm = tm;
            this.region = region;
        }
    }

    /**
     * @param genomeKmerSet set of encoded k-mer hashes from the reference
     * genome
     * @param sequence target DNA sequence
     * @param sequenceName human-readable name/ID for the target
     * @param minPrimerLength minimum acceptable primer length (nt)
     * @param maxPrimerLength maximum acceptable primer length (nt)
     * @param minMeltingTemp minimum melting temperature (°C)
     * @param maxMeltingTemp maximum melting temperature (°C)
     * @param minLinguisticComplexity minimum linguistic complexity score
     * @param exonBoundaries pairs of [start, end) positions for each exon
     * @param existingPrimers previously designed primers to check for
     * cross-dimers
     * @param referenceSeq reference sequence used for masking
     * @param fastaSeq additional FASTA sequence for homology masking
     * @param useHomologyFilter whether to apply separate homology-based
     * filtering
     */
    public PrimerDesign(LongIntHashMap genomeKmerSet, String sequence, String sequenceName,
            int minPrimerLength, int maxPrimerLength,
            int minMeltingTemp, int maxMeltingTemp,
            int minLinguisticComplexity, int[] exonBoundaries,
            String[] existingPrimers, String referenceSeq,
            String fastaSeq, boolean useHomologyFilter) {
        this.genomeKmerSet = genomeKmerSet;
        this.sequence = sequence;
        this.sequenceName = sequenceName;
        this.referenceSeq = referenceSeq;
        this.fastaSeq = fastaSeq;
        this.minLinguisticComplexity = minLinguisticComplexity;
        this.minPrimerLength = minPrimerLength;
        this.maxPrimerLength = maxPrimerLength;
        this.minMeltingTemp = minMeltingTemp;
        this.maxMeltingTemp = maxMeltingTemp;
        this.exonBoundaries = exonBoundaries;
        this.existingPrimers = existingPrimers;
        this.useHomologyFilter = useHomologyFilter;
    }

    public PrimersCollector[] getForwardPrimers() {
        return forwardPrimers;
    }

    public PrimersCollector[] getReversePrimers() {
        return reversePrimers;
    }

    /**
     * Per region, the widened search flank (bp) with which the region was
     * rescued, or 0 where the default flank was used.
     */
    public int[] getWidenedFlanks() {
        return widenedFlanks;
    }

    /**
     * Every valid combination of the designed primers, a superset of the pairs
     * returned by {@link #runDesign}; the panel falls back on them to close gaps.
     */
    public List<Pair> getAllPairs() {
        return allPairs;
    }

    /** What the primer search found, for the panel map; set by {@link #runDesign}. */
    public SearchMap getSearchMap() {
        return search;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Genome masking
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * Marks positions in {@code mask} that are covered by k-mers present in the
     * genome reference set, scanning both the forward and reverse-complement
     * strands.
     */
    private int[] applyGenomeMasking(int[] mask, String seq, int kmerSize) {
        int seqLen = seq.length();
        byte[] fwdBases = encodeBases(seq);
        maskStrand(mask, fwdBases, seqLen, kmerSize, /* reverseMapping */ false);

        fwdBases = encodeBases(dna.ComplementDNA(seq));
        maskStrand(mask, fwdBases, seqLen, kmerSize, /* reverseMapping */ true);

        return mask;
    }

    /**
     * Converts a DNA string to its numeric base encoding via
     * {@code tables.dx2}.
     */
    private byte[] encodeBases(String dnaStr) {
        byte[] bases = dnaStr.getBytes();
        for (int i = 0; i < bases.length; i++) {
            bases[i] = tables.dx2[bases[i]];
        }
        return bases;
    }

    /**
     * Slides a k-mer window across one strand, incrementing mask positions
     * wherever the k-mer is found in {@link #genomeKmerSet}.
     *
     * @param reverseMapping if true, map hits back to the forward-strand
     * coordinates
     */
    private void maskStrand(int[] mask, byte[] bases, int seqLen, int kmerSize, boolean reverseMapping) {
        int[] window = new int[kmerSize];
        int[] baseCounts = new int[5]; // A, C, G, T, N (index 4 = ambiguous)

        // Seed the sliding window with the first (kmerSize-1) bases
        for (int i = 0; i < kmerSize - 1; i++) {
            window[i] = bases[i];
            baseCounts[window[i]]++;
        }

        for (int i = kmerSize - 1; i < seqLen; i++) {
            window[kmerSize - 1] = bases[i];
            baseCounts[window[kmerSize - 1]]++;

            int windowStart = i - kmerSize + 1;
            boolean noAmbiguousBases = (baseCounts[4] == 0);

            if (noAmbiguousBases) {
                long kmerCode = encodeKmer(bases, windowStart, kmerSize);
                if (genomeKmerSet.get(kmerCode) > 1) {
                    int maskStart = reverseMapping ? (seqLen - windowStart - kmerSize) : windowStart;
                    int maskEnd = maskStart + kmerSize;
                    for (int j = maskStart; j < maskEnd; j++) {
                        mask[j]++;
                    }
                }
            }

            // Slide: remove the leftmost base, shift window
            baseCounts[bases[i + 1 - kmerSize]]--;
            System.arraycopy(window, 1, window, 0, kmerSize - 1);
        }
    }

    /**
     * Encodes a k-mer as a long by treating each numeric base as a decimal
     * digit. For k ≤ 19 the result is always non-negative; larger k values may
     * overflow into negative territory but remain unique.
     */
    private long encodeKmer(byte[] bases, int start, int kmerSize) {
        long code = 0L;
        for (int i = start; i < start + kmerSize; i++) {
            code = code * 5 + bases[i];
        }
        return code;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Main design pipeline
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * Runs the full primer-design pipeline: mask the sequence, merge exons,
     * design forward and reverse primers for each exon region, and combine them
     * into valid PCR pairs.
     *
     * @param forwardTail 5′ tail appended to forward primers
     * @param reverseTail 5′ tail appended to reverse primers
     * @param required5PrimeEnd required pattern at the 5′ end of primers
     * @param required3PrimeEnd required pattern at the 3′ end of primers
     * @param primerOverlap allowed overlap between adjacent primer tiles (nt)
     * @param minPcrProduct minimum acceptable PCR product size (bp)
     * @param maxPcrProduct maximum acceptable PCR product size (bp)
     * @return list of valid forward–reverse primer pairs
     */
    public List<Pair> runDesign(String forwardTail, String reverseTail, String required5PrimeEnd, String required3PrimeEnd, int primerOverlap, int minPcrProduct, int maxPcrProduct) {

        int seqLen = sequence.length();

        // ── Build masking arrays ─────────────────────────────────────────
        int[] senseMask;
        int[] homologyMask = new int[seqLen];

        if (useHomologyFilter) {
            senseMask = new MaskingSequence(GENOME_KMER_SIZE).mask(sequence + referenceSeq);
            homologyMask = new MaskingSequence(GENOME_KMER_SIZE).mask(sequence + fastaSeq);
        } else {
            senseMask = new MaskingSequence(GENOME_KMER_SIZE).mask(sequence + referenceSeq + fastaSeq);
        }

        senseMask = Arrays.copyOf(senseMask, seqLen);
        int[] genomeMask = applyGenomeMasking(new int[seqLen], sequence, GENOME_KMER_SIZE);
        search = new SearchMap(genomeKmerSet != null && !genomeKmerSet.isEmpty());
        search.repeats.addAll(stretches(senseMask));
        search.genome.addAll(stretches(genomeMask));
        for (int i = 0; i < seqLen; i++) {
            senseMask[i] += genomeMask[i];
        }

        // Antisense masks are the sense masks in reverse order
        int[] antisenseMask = reverseArray(senseMask, seqLen);
        int[] homologyMaskReverse = reverseArray(homologyMask, seqLen);

        // ── Merge overlapping / close exons ──────────────────────────────
        int[] targetBoundaries = exonBoundaries.clone(); // the regions to cover, before merging
        search.lowComplexity.addAll(lowComplexity(targetBoundaries, maxPcrProduct));
        mergeExons(minPcrProduct, maxPcrProduct);

        if (exonBoundaries.length < 2) {
            return new ArrayList<>();
        }

        // ── Design primers per exon region ───────────────────────────────
        int extraSeq = Math.min(minPcrProduct, MAX_EXTRA_SEQ);
        int numExons = exonBoundaries.length / 2;
        forwardPrimers = new PrimersCollector[numExons];
        reversePrimers = new PrimersCollector[numExons];
        widenedFlanks = new int[numExons];
        PCRcollector pcrCollector = new PCRcollector();
        String antisenseSeq = dna.ComplementDNA(sequence);

        int taskIndex = 0;
        for (int i = 0; i < exonBoundaries.length - 1; i += 2) {
            int exonStart = exonBoundaries[i];
            int exonEnd = exonBoundaries[i + 1];
            int exonLength = exonEnd - exonStart;

            // Merged into a preceding region, whose primers cover it
            if (exonStart < 0) {
                forwardPrimers[taskIndex] = buildCollector(new ArrayList<>(), true, seqLen);
                reversePrimers[taskIndex] = buildCollector(new ArrayList<>(), false, seqLen);
                setTaskMetadata(forwardPrimers[taskIndex], taskIndex + 1, exonStart, exonEnd);
                setTaskMetadata(reversePrimers[taskIndex], taskIndex + 1, exonStart, exonEnd);
                search.pairing.add(null);
                taskIndex++;
                continue;
            }

            designRegionPrimers(taskIndex, exonStart, exonEnd, extraSeq, extraSeq, antisenseSeq, senseMask, homologyMask, antisenseMask, homologyMaskReverse, forwardTail, reverseTail, required5PrimeEnd, required3PrimeEnd, primerOverlap);

            // ── Rescue ───────────────────────────────────────────────────
            // When the primers found in the default flanks cannot cover the region
            // (no pair at all, or a region edge out of reach), search again as far
            // from the region as maxPCR allows, so that it is not dropped from the panel.
            List<int[]> targets = regionTargets(taskIndex, targetBoundaries);
            int covered = candidateCoverage(taskIndex, targets, minPcrProduct, maxPcrProduct);
            int wideFlank = maxPcrProduct - Math.min(exonLength, extraSeq) - minPrimerLength;
            boolean roomToWiden = exonStart - extraSeq > 0 || exonEnd + extraSeq < seqLen;
            if (covered < targetLength(targets) && wideFlank > extraSeq && roomToWiden) {
                PrimersCollector defaultForward = forwardPrimers[taskIndex];
                PrimersCollector defaultReverse = reversePrimers[taskIndex];
                designRegionPrimers(taskIndex, exonStart, exonEnd, extraSeq, wideFlank, antisenseSeq, senseMask, homologyMask, antisenseMask, homologyMaskReverse, forwardTail, reverseTail, required5PrimeEnd, required3PrimeEnd, primerOverlap);
                if (candidateCoverage(taskIndex, targets, minPcrProduct, maxPcrProduct) > covered) {
                    widenedFlanks[taskIndex] = wideFlank;
                } else {
                    forwardPrimers[taskIndex] = defaultForward;
                    reversePrimers[taskIndex] = defaultReverse;
                }
            }

            addSites(taskIndex);
            search.pairing.add(pairingNote(taskIndex, exonStart, exonEnd, minPcrProduct, maxPcrProduct));

            // ── Collect valid pair ───────────────────────────────────────
            if (forwardPrimers[taskIndex].Amount() > 0 && reversePrimers[taskIndex].Amount() > 0) {
                pcrCollector.add(forwardPrimers[taskIndex], reversePrimers[taskIndex]);
            }

            taskIndex++;
        }

        if (pcrCollector.getCurrentSize() == 0) {
            return new ArrayList<>();
        }
        allPairs = pcrCollector.CombineAllPrimers(minPcrProduct, maxPcrProduct);
        return pcrCollector.CombinePrimers(minPcrProduct, maxPcrProduct);
    }

    /**
     * Designs the forward and reverse primers of one region, searching up to
     * {@code flank} bp outside it. A region shorter than {@code extraSeq} needs an
     * amplicon that spans it, so its primers are searched outside it only.
     */
    private void designRegionPrimers(int taskIndex, int exonStart, int exonEnd, int extraSeq, int flank, String antisenseSeq, int[] senseMask, int[] homologyMask, int[] antisenseMask, int[] homologyMaskReverse, String forwardTail, String reverseTail, String required5PrimeEnd, String required3PrimeEnd, int primerOverlap) {
        int seqLen = sequence.length();
        int exonLength = exonEnd - exonStart;

        // ── Forward primers ──────────────────────────────────────────
        int fwdSearchStart = Math.max(0, exonStart - flank);
        int fwdSearchEnd = (exonLength < extraSeq) ? exonStart : exonEnd;
        fwdSearchEnd = Math.min(fwdSearchEnd, seqLen - 1);

        forwardPrimers[taskIndex] = designPrimersForRegion(true, sequence, forwardTail, senseMask, homologyMask, minLinguisticComplexity, minPrimerLength, maxPrimerLength, minMeltingTemp, maxMeltingTemp, fwdSearchStart, fwdSearchEnd, primerOverlap, sequenceName, required5PrimeEnd, required3PrimeEnd);

        if (forwardPrimers[taskIndex].Amount() < MIN_PRIMERS_BEFORE_RELAXING) {
            forwardPrimers[taskIndex] = designPrimersForRegion(true, sequence, forwardTail, senseMask, homologyMask, minLinguisticComplexity - RELAXED_LC_OFFSET, minPrimerLength, maxPrimerLength + RELAXED_MAXLEN_OFFSET, minMeltingTemp, maxMeltingTemp, fwdSearchStart, fwdSearchEnd, primerOverlap, sequenceName, "", "");
        }
        setTaskMetadata(forwardPrimers[taskIndex], taskIndex + 1, exonStart, exonEnd);

        // ── Reverse primers ──────────────────────────────────────────
        int revSearchStart = seqLen - exonEnd - flank;
        int revSearchEnd = (exonLength < extraSeq) ? seqLen - exonEnd : seqLen - exonStart + flank;
        revSearchStart = Math.max(0, revSearchStart);
        revSearchEnd = Math.min(revSearchEnd, seqLen - 1);

        reversePrimers[taskIndex] = designPrimersForRegion(false, antisenseSeq, reverseTail, antisenseMask, homologyMaskReverse, minLinguisticComplexity, minPrimerLength, maxPrimerLength, minMeltingTemp, maxMeltingTemp, revSearchStart, revSearchEnd, primerOverlap, sequenceName, required5PrimeEnd, required3PrimeEnd);

        if (reversePrimers[taskIndex].Amount() < MIN_PRIMERS_BEFORE_RELAXING) {
            reversePrimers[taskIndex] = designPrimersForRegion(false, antisenseSeq, reverseTail, antisenseMask, homologyMaskReverse, minLinguisticComplexity - RELAXED_LC_OFFSET, minPrimerLength, maxPrimerLength + RELAXED_MAXLEN_OFFSET, minMeltingTemp, maxMeltingTemp, revSearchStart, revSearchEnd, primerOverlap, sequenceName, "", "");
        }
        setTaskMetadata(reversePrimers[taskIndex], taskIndex + 1, exonStart, exonEnd);
    }

    /**
     * The original (pre-merge) regions that region {@code taskIndex} has to cover,
     * in the amplicon frame: the region itself and every region merged into it.
     */
    private List<int[]> regionTargets(int taskIndex, int[] targetBoundaries) {
        List<int[]> targets = new ArrayList<>();
        for (int k = 0; k + 1 < exonBoundaries.length; k += 2) {
            if (k == 2 * taskIndex || exonBoundaries[k] == -(taskIndex + 1)) {
                int[] t = CoverageReport.targetInterval(targetBoundaries[k], targetBoundaries[k + 1]);
                if (t != null) {
                    targets.add(t);
                }
            }
        }
        return targets;
    }

    private static int targetLength(List<int[]> targets) {
        int n = 0;
        for (int[] t : targets) {
            n += t[1] - t[0];
        }
        return n;
    }

    /** Records the candidate primers of region {@code taskIndex} for the panel map. */
    private void addSites(int taskIndex) {
        PrimersCollector fp = forwardPrimers[taskIndex];
        for (int i = 0; i < fp.Amount(); i++) {
            int x = fp.getPrimerLocations()[i];
            search.sites.add(new Site(true, x, x + fp.getPrimerLengths()[i], fp.getpPrimerName()[i], fp.getTms()[i], taskIndex));
        }
        PrimersCollector rp = reversePrimers[taskIndex];
        for (int i = 0; i < rp.Amount(); i++) {
            int x = rp.getPrimerLocations()[i]; // 5' end, the right end of the binding site
            search.sites.add(new Site(false, x - rp.getPrimerLengths()[i] + 1, x + 1, rp.getpPrimerName()[i], rp.getTms()[i], taskIndex));
        }
    }

    /**
     * Why region {@code taskIndex} gets no primer pair, or null when it gets
     * one: the pairing rules of {@link PCRcollector#pairAllPrimers} applied to
     * its candidate primers.
     */
    private String pairingNote(int taskIndex, int exonStart, int exonEnd, int minPcr, int maxPcr) {
        PrimersCollector fp = forwardPrimers[taskIndex];
        PrimersCollector rp = reversePrimers[taskIndex];
        int nf = fp.Amount();
        int nr = rp.Amount();
        if (nf == 0 && nr == 0) {
            return "no forward or reverse primer passes the filters";
        }
        if (nf == 0) {
            return "no forward primer passes the filters (" + nr + " reverse)";
        }
        if (nr == 0) {
            return "no reverse primer passes the filters (" + nf + " forward)";
        }
        int exon1 = exonStart;
        int exon2 = exonEnd;
        if (exon2 - exon1 < minPcr) { // a short region needs an amplicon spanning it
            exon1 = exonEnd;
            exon2 = exonStart;
        }
        int[] fx = fp.getPrimerLocations();
        int[] rx = rp.getPrimerLocations();
        int longestBelow = -1;              // amplicon sizes as reported, 1 + rx - fx
        int shortestAbove = Integer.MAX_VALUE;
        int inRange = 0;
        int dimers = 0;
        for (int f = 0; f < nf; f++) {
            for (int r = 0; r < nr; r++) {
                int tpcr = rx[r] - fx[f];
                if (rx[r] <= 0 || rx[r] <= exon1 || fx[f] >= exon2 || tpcr <= 0) {
                    continue;
                }
                if (tpcr < minPcr) {
                    longestBelow = Math.max(longestBelow, tpcr + 1);
                } else if (tpcr > maxPcr) {
                    shortestAbove = Math.min(shortestAbove, tpcr + 1);
                } else {
                    inRange++;
                    if (Oligo.quickDimer(fp.getPrimer()[f], rp.getPrimer()[r]) != 0) {
                        dimers++;
                    }
                }
            }
        }
        if (inRange > dimers) {
            return null;
        }
        if (inRange > 0) {
            return "all " + inRange + " pairs of " + minPcr + "-" + maxPcr + " bp form a primer-dimer";
        }
        String found = nf + " forward and " + nr + " reverse primers, but ";
        if (longestBelow > 0 && shortestAbove < Integer.MAX_VALUE) {
            return found + "their pairs give " + longestBelow + " bp or less, or " + shortestAbove + " bp or more, none within minPCR-maxPCR (" + minPcr + "-" + maxPcr + " bp)";
        }
        if (shortestAbove < Integer.MAX_VALUE) {
            return found + "the shortest pair gives " + shortestAbove + " bp, longer than maxPCR (" + maxPcr + " bp)";
        }
        if (longestBelow > 0) {
            return found + "the longest pair gives " + longestBelow + " bp, shorter than minPCR (" + minPcr + " bp)";
        }
        return found + "no forward primer lies upstream of a reverse primer around the region";
    }

    /**
     * Target bases of the region that its primer pairs can cover between them.
     */
    private int candidateCoverage(int taskIndex, List<int[]> targets, int minPcrProduct, int maxPcrProduct) {
        List<int[]> amplicons = new ArrayList<>();
        for (Pair p : PCRcollector.pairAllPrimers(forwardPrimers[taskIndex], reversePrimers[taskIndex], minPcrProduct, maxPcrProduct)) {
            amplicons.add(CoverageReport.interval(p.fx5, 1 + p.rx5));
        }
        List<int[]> merged = CoverageReport.merge(amplicons);
        int n = 0;
        for (int[] t : targets) {
            n += CoverageReport.overlap(merged, t[0], t[1]);
        }
        return n;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Primer candidate evaluation
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * Candidate primer, collecting all scored properties in one place instead
     * of parallel arrays.
     */
    private static class PrimerCandidate {

        final String fullSequence; // tail + binding region
        final String label;       // human-readable name
        final int position;
        final int bindingLength;
        final double meltingTemp;
        final double gcContent;
        final int complexity;

        PrimerCandidate(String fullSequence, String label, int position, int bindingLength, double meltingTemp, double gcContent, int complexity) {
            this.fullSequence = fullSequence;
            this.label = label;
            this.position = position;
            this.bindingLength = bindingLength;
            this.meltingTemp = meltingTemp;
            this.gcContent = gcContent;
            this.complexity = complexity;
        }
    }

    /**
     * Designs primer candidates for a single region of the (sense or antisense)
     * strand. Walks a sliding window, keeping a running count of masked and
     * homology-flagged positions so that only "clean" windows are evaluated.
     */
    private PrimersCollector designPrimersForRegion(boolean isForward, String strandSeq, String tail, int[] mask, int[] homologyMask, int minLC, int minLen, int maxLen, double minTm, double maxTm, int searchStart, int searchEnd, int primerOverlap, String name, String end5Pattern, String end3Pattern) {

        int seqLen = strandSeq.length();
        List<PrimerCandidate> candidates = new ArrayList<>();

        VariantGenerator gen3 = new VariantGenerator(end3Pattern);
        String allowed3 = gen3.generatedstr();
        int req3Len = gen3.minstrlen();

        VariantGenerator gen5 = new VariantGenerator(end5Pattern);
        String allowed5 = gen5.generatedstr();
        int req5Len = gen5.minstrlen();

        primer calculator = new primer();

        // Running counts inside the current minimum-length window
        int maskedCount = 0;
        int homologyCount = 0;

        // Seed window [searchStart .. searchStart+minLen-2]
        for (int y = searchStart; y < searchStart + minLen - 1; y++) {
            if (mask[y] > 0) {
                maskedCount++;
            }
            if (useHomologyFilter && homologyMask[y] == 0) {
                homologyCount++;
            }
        }

        // Slide across the search region
        for (int pos = searchStart; pos < searchEnd - minLen + 1; pos++) {
            // Extend window by one base on the right
            if (mask[pos + minLen - 1] > 0) {
                maskedCount++;
            }
            if (useHomologyFilter && homologyMask[pos + minLen - 1] == 0) {
                homologyCount++;
            }

            if (maskedCount <= MASK_TOLERANCE && homologyCount == 0) {
                PrimerCandidate candidate = evaluateCandidate(strandSeq, tail, pos, minLen, maxLen, searchStart, minTm, maxTm, minLC, allowed3, req3Len, allowed5, req5Len, calculator, isForward, seqLen, name);

                if (candidate != null) {
                    addCandidateIfBetter(candidates, candidate, primerOverlap);
                }
            }

            // Shrink window by removing the leftmost base
            if (mask[pos] > 0) {
                maskedCount--;
            }
            if (useHomologyFilter && homologyMask[pos] == 0) {
                homologyCount--;
            }
        }

        return buildCollector(candidates, isForward, seqLen);
    }

    /**
     * Evaluates a single candidate position: extends the primer to reach the
     * target Tm, checks 3′/5′-end constraints, SSRs, complexity, and dimers.
     *
     * @return a scored {@link PrimerCandidate}, or {@code null} if the position
     * fails any filter
     */
    private PrimerCandidate evaluateCandidate(String strandSeq, String tail, int anchorPos, int minLen, int maxLen, int searchStart, double minTm, double maxTm, int minLC, String allowed3, int req3Len, String allowed5, int req5Len, primer calculator, boolean isForward, int seqLen, String name) {

        // Check 3′-end constraint
        String threePrimeEnd = strandSeq.substring(anchorPos + minLen - req3Len, anchorPos + minLen);
        if (req3Len > 0 && !allowed3.contains(threePrimeEnd)) {
            return null;
        }

        // Extend 5′ until we reach the minimum Tm (or run out of room)
        double tm = 0;
        String bindingRegion = null;
        int primerStart = anchorPos;

        for (int extension = 0; extension <= maxLen - minLen; extension++) {
            primerStart = anchorPos - extension;
            if (primerStart < searchStart) {
                break;
            }

            bindingRegion = strandSeq.substring(primerStart, primerStart + minLen + extension);
            if (allowed5.contains(bindingRegion.substring(0, req5Len))) {
                tm = calculator.getTm(bindingRegion);
            }
            if (tm >= minTm) {
                break;
            }
        }

        if (tm < minTm || tm > maxTm) {
            return null;
        }
        if (bindingRegion == null || containsSSR(bindingRegion)) {
            return null;
        }

        int complexity = calculator.getLC();
        if (complexity < minLC) {
            return null;
        }

        // Dimer checks
        String fullPrimer = tail + bindingRegion;
        if (Oligo.quickDimer(fullPrimer) != 0) {
            return null;
        }
        if (existingPrimers.length > 0 && Oligo.DimersCheck(bindingRegion, existingPrimers) != 0) {
            return null;
        }

        // Build label and position
        int reportedPos;
        String label;
        if (isForward) {
            reportedPos = primerStart;
            label = name + ":F_" + (primerStart + 1) + "-" + (primerStart + bindingRegion.length());
        } else {
            reportedPos = primerStart; // will be flipped later in buildCollector
            int absPos = seqLen - primerStart - 1;
            label = name + ":R_" + (absPos + 1) + "-" + (seqLen - primerStart - bindingRegion.length() + 1);
        }

        return new PrimerCandidate(
                fullPrimer, label, reportedPos,
                bindingRegion.length(), tm, calculator.getGC(), complexity);
    }

    /**
     * Adds a candidate to the list, enforcing the overlap/tiling constraint: if
     * the new candidate overlaps the last one, keep whichever has higher
     * linguistic complexity.
     */
    private void addCandidateIfBetter(List<PrimerCandidate> candidates, PrimerCandidate candidate, int primerOverlap) {
        if (!candidates.isEmpty()) {
            PrimerCandidate last = candidates.get(candidates.size() - 1);
            int nextAllowedPos = last.position + last.bindingLength - primerOverlap;

            if (nextAllowedPos >= candidate.position) {
                // Overlaps — keep the one with better complexity
                if (last.complexity < candidate.complexity) {
                    candidates.remove(candidates.size() - 1);
                } else {
                    return; // existing one is better or equal
                }
            }
        }
        candidates.add(candidate);
    }

    /**
     * Converts the candidate list into a {@link PrimersCollector}, flipping
     * positions for reverse-strand primers.
     */
    private PrimersCollector buildCollector(List<PrimerCandidate> candidates, boolean isForward, int seqLen) {
        int n = candidates.size();
        String[] sequences = new String[n];
        String[] names = new String[n];
        int[] lengths = new int[n];
        double[] temps = new double[n];
        double[] gcContents = new double[n];
        int[] complexities = new int[n];
        int[] positions = new int[n];

        for (int i = 0; i < n; i++) {
            PrimerCandidate c = candidates.get(i);
            sequences[i] = c.fullSequence;
            names[i] = c.label;
            lengths[i] = c.bindingLength;
            temps[i] = c.meltingTemp;
            gcContents[i] = c.gcContent;
            complexities[i] = c.complexity;
            positions[i] = isForward ? c.position : (seqLen - c.position - 1);
        }

        PrimersCollector pc = new PrimersCollector();
        pc.insert(n, sequences, names, lengths, temps, gcContents, complexities, positions);
        return pc;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Helpers
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * Merges exon pairs that can be covered by a single PCR product or are too
     * close together to amplify independently. Modifies {@link #exonBoundaries}
     * in place (marks merged entries with 0/-index).
     */
    private void mergeExons(int minPcrProduct, int maxPcrProduct) {
        for (int i = 0; i < exonBoundaries.length - 1; i += 2) {
            if (exonBoundaries[i + 1] <= 0) {
                continue;
            }

            for (int j = i + 2; j < exonBoundaries.length - 1; j += 2) {
                if (exonBoundaries[j + 1] <= 0) {
                    continue;
                }

                boolean fitsInOnePcr = (exonBoundaries[j + 1] - exonBoundaries[i]) < maxPcrProduct;
                boolean tooCloseToSplit = (exonBoundaries[j] - exonBoundaries[i + 1]) < minPcrProduct;

                if (fitsInOnePcr || tooCloseToSplit) {
                    // Absorb exon j into exon i; j may end inside i (overlapping annotations)
                    exonBoundaries[i + 1] = Math.max(exonBoundaries[i + 1], exonBoundaries[j + 1]);
                    exonBoundaries[j + 1] = 0;
                    exonBoundaries[j] = -((i + 2) / 2);
                } else {
                    break;
                }
            }
        }
    }

    /**
     * Returns a new array with elements in reverse order.
     */
    private int[] reverseArray(int[] src, int length) {
        int[] reversed = Arrays.copyOf(src, length);
        for (int i = 0, j = length - 1; i < j; i++, j--) {
            int tmp = reversed[i];
            reversed[i] = reversed[j];
            reversed[j] = tmp;
        }
        return reversed;
    }

    private void setTaskMetadata(PrimersCollector pc, int taskId, int x1, int x2) {
        pc.setTaskId(taskId);
        pc.setTaskX1(x1);
        pc.setTaskX2(x2);
    }

    /**
     * Checks whether the primer contains any disqualifying simple-sequence
     * repeats.
     */
    private static boolean containsSSR(String primer) {
        for (String pattern : SSR_PATTERNS) {
            if (primer.contains(pattern)) {
                return true;
            }
        }
        return false;
    }

    /** The stretches where {@code mask} is set, 0-based half-open. */
    private static List<int[]> stretches(int[] mask) {
        List<int[]> out = new ArrayList<>();
        int start = -1;
        for (int i = 0; i <= mask.length; i++) {
            boolean set = i < mask.length && mask[i] > 0;
            if (set && start < 0) {
                start = i;
            } else if (!set && start >= 0) {
                out.add(new int[]{start, i});
                start = -1;
            }
        }
        return out;
    }

    /**
     * Bases within {@code reach} bp of the targets that no primer of the
     * minimum length can include, because every such window falls below the
     * minimum linguistic complexity or holds a simple sequence repeat; 0-based
     * half-open stretches. Only the surroundings of the targets are checked, as
     * only there can their primers lie.
     */
    private List<int[]> lowComplexity(int[] targetBoundaries, int reach) {
        int seqLen = sequence.length();
        List<int[]> areas = new ArrayList<>();
        for (int k = 0; k + 1 < targetBoundaries.length; k += 2) {
            int[] t = CoverageReport.targetInterval(targetBoundaries[k], targetBoundaries[k + 1]);
            if (t != null) {
                areas.add(new int[]{Math.max(0, t[0] - reach), Math.min(seqLen, t[1] + reach)});
            }
        }
        primer calculator = new primer();
        List<int[]> out = new ArrayList<>();
        for (int[] a : CoverageReport.merge(areas)) {
            // windows reach past the area, so that its edge bases are judged like the others
            int from = Math.max(0, a[0] - minPrimerLength);
            int to = Math.min(seqLen, a[1] + minPrimerLength);
            int windows = to - from - minPrimerLength + 1;
            if (windows <= 0) {
                continue;
            }
            int[] passing = new int[windows + 1]; // passing[w]: windows before w that pass
            for (int w = 0; w < windows; w++) {
                String s = sequence.substring(from + w, from + w + minPrimerLength);
                boolean pass = calculator.getLC(s) >= minLinguisticComplexity && !containsSSR(s);
                passing[w + 1] = passing[w] + (pass ? 1 : 0);
            }
            int start = -1;
            for (int i = a[0]; i <= a[1]; i++) {
                boolean low = false;
                if (i < a[1]) {
                    int lo = Math.max(0, i - from - minPrimerLength + 1); // windows holding base i
                    int hi = Math.min(windows - 1, i - from);
                    low = hi < lo || passing[hi + 1] - passing[lo] == 0;
                }
                if (low && start < 0) {
                    start = i;
                } else if (!low && start >= 0) {
                    out.add(new int[]{start, i});
                    start = -1;
                }
            }
        }
        return out;
    }
}
