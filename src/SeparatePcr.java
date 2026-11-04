import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Primer pairs for what the panel leaves uncovered, to be amplified in a PCR of
 * their own.
 *
 * A pair in its own tube has to get along with nothing but itself, so what
 * comes from sharing a tube falls away: its Tm need only match that of its
 * partner, and no other primer of the panel is checked. Each uncovered stretch
 * is searched with the filters of the panel first and then, step by step, with
 * looser ones; the first step that gives a pair is used, and the report says
 * what it relaxed.
 *
 * The primers lie outside what they have to amplify, by a margin, so that the
 * uncovered bases are read clear of the primers' own sequence. A stretch too
 * long for one amplicon is tiled, and neighbouring stretches share an amplicon
 * where they can.
 *
 * Coordinates are 0-based half-open on the sense strand.
 */
public class SeparatePcr {

    private static final int PAIRS = 3;        // pairs offered per amplicon: the best one and its alternatives
    private static final int MIN_SHIFT = 10;   // bp between the binding sites of the pairs offered
    private static final int SKIP = 10;        // bp to move along when no pair reaches a base
    private static final int MARGIN = 30;      // bp between a primer and what it has to amplify
    private static final int LAST_MARGIN = 10; // the same, on the last step
    private static final int CLEAN_3 = 8;      // bases at the 3' end that may not be masked
    private static final int UNIQUE_3 = 15;    // bases at the 3' end that have to be unique in the target sequence
    private static final int MIN_PRIMER = 15;
    private static final int MIN_AMPLICON = 100;
    private static final int LAST_MASKED = 5;  // masked bases tolerated on the last step
    private static final int LINE = 60;

    /** The design settings of the run. */
    public static final class Settings {

        final int minPcr;
        final int maxPcr;
        final int minLen;
        final int maxLen;
        final int minTm;
        final int maxTm;
        final int minLc;
        final int maxSeparate;
        final String e5;
        final String e3;
        final String fTail;
        final String rTail;

        /**
         * @param separatePcr largest amplicon of a separate PCR; 0 for one and a
         *                    half times maxPcr
         */
        public Settings(int minPcr, int maxPcr, int minLen, int maxLen, int minTm, int maxTm, int minLc, int separatePcr,
                String e5, String e3, String fTail, String rTail) {
            this.minPcr = minPcr;
            this.maxPcr = maxPcr;
            this.minLen = minLen;
            this.maxLen = maxLen;
            this.minTm = minTm;
            this.maxTm = maxTm;
            this.minLc = minLc;
            this.maxSeparate = (separatePcr > 0) ? Math.max(separatePcr, maxPcr) : maxPcr * 3 / 2;
            this.e5 = e5;
            this.e3 = e3;
            this.fTail = fTail;
            this.rTail = rTail;
        }
    }

    /** One step of the search: the filters a primer pair has to pass. */
    private static final class Step {

        String relaxed = "nothing, these are the filters of the panel";
        int minTm;
        int maxTm;
        double maxTmDifference;
        int minLen;
        int maxLen;
        boolean anyEnds = false;
        int minLc;
        int maxMasked = PrimerDesign.MASK_TOLERANCE;
        int minSize;
        int maxSize;
        int margin = MARGIN;
        final Map<Integer, Primer> forward = new HashMap<>(); // by 3' end: primer, or null where none passes
        final Map<Integer, Primer> reverse = new HashMap<>();
    }

    /** A primer that passed the filters of a step. */
    private static final class Primer {

        int start;
        int end;
        String name;
        String binding;
        String full;        // with the adapter tail
        double tm;
        double gc;
        int lc;
        int masked;
        int hits3 = -1;     // occurrences of its 3' end in the target sequence; -1 = not counted yet
    }

    /** A forward and a reverse primer, and the bases they amplify with the margin kept. */
    private static final class Choice {

        Primer f;
        Primer r;
        int lo;
        int hi;
        int gain;           // uncovered bases it amplifies

        int size() {
            return r.end - f.start;
        }

        double tmDifference() {
            return Math.abs(f.tm - r.tm);
        }
    }

    /** One amplicon: the best pair first, then alternatives that amplify the same bases. */
    private static final class Amplicon {

        Step step;
        List<Choice> pairs;
        int lo;
        int hi;
    }

    private final String name;
    private final String sequence;
    private final String antisense;
    private final int seqLen;
    private final boolean[] masked;
    private final boolean genomeChecked;
    private final Settings c;
    private final List<Step> steps;
    private final String allowed3;
    private final int req3;
    private final String allowed5;
    private final int req5;
    private final primer calculator = new primer();

    /**
     * @param name     record ID, as in the primer names
     * @param sequence the target sequence
     * @param search   where primers cannot bind, as found by the primer design
     */
    public SeparatePcr(String name, String sequence, PrimerDesign.SearchMap search, Settings settings) {
        this.name = name;
        this.sequence = sequence;
        this.antisense = dna.ComplementDNA(sequence);
        this.seqLen = sequence.length();
        this.c = settings;
        this.masked = new boolean[seqLen];
        mark(search.repeats);
        mark(search.genome);
        this.genomeChecked = search.genomeChecked;
        this.steps = steps(settings);
        VariantGenerator g3 = new VariantGenerator(settings.e3);
        this.allowed3 = g3.generatedstr();
        this.req3 = g3.minstrlen();
        VariantGenerator g5 = new VariantGenerator(settings.e5);
        this.allowed5 = g5.generatedstr();
        this.req5 = g5.minstrlen();
    }

    private void mark(List<int[]> stretches) {
        for (int[] v : stretches) {
            Arrays.fill(masked, v[0], v[1], true);
        }
    }

    /** The filters of the panel, then looser ones. Each step keeps what the steps before it relaxed. */
    private static List<Step> steps(Settings c) {
        int tm3Low = c.minTm - 3;
        int tm3High = c.maxTm + 3;
        int shortest = Math.min(c.minLen, Math.max(MIN_PRIMER, c.minLen - 2));
        int longest = c.maxLen + 4;
        int smallest = Math.min(c.minPcr, MIN_AMPLICON);

        Step panel = new Step();
        panel.minTm = c.minTm;
        panel.maxTm = c.maxTm;
        panel.maxTmDifference = c.maxTm - c.minTm;
        panel.minLen = c.minLen;
        panel.maxLen = c.maxLen;
        panel.minLc = c.minLc;
        panel.minSize = c.minPcr;
        panel.maxSize = c.maxPcr;

        Step tm = new Step();
        tm.relaxed = "Tm " + tm3Low + "-" + tm3High + " °C";
        tm.minTm = tm3Low;
        tm.maxTm = tm3High;
        tm.maxTmDifference = 3;
        tm.minLen = c.minLen;
        tm.maxLen = c.maxLen;
        tm.minLc = c.minLc;
        tm.minSize = c.minPcr;
        tm.maxSize = c.maxPcr;

        Step ends = new Step();
        ends.relaxed = tm.relaxed + ", primers " + shortest + "-" + longest + " nt, any 3' and 5' end";
        ends.minTm = tm3Low;
        ends.maxTm = tm3High;
        ends.maxTmDifference = 3;
        ends.minLen = shortest;
        ends.maxLen = longest;
        ends.anyEnds = true;
        ends.minLc = c.minLc;
        ends.minSize = c.minPcr;
        ends.maxSize = c.maxPcr;

        Step size = new Step();
        size.relaxed = ends.relaxed + ", amplicons " + smallest + "-" + c.maxSeparate + " bp";
        size.minTm = tm3Low;
        size.maxTm = tm3High;
        size.maxTmDifference = 3;
        size.minLen = shortest;
        size.maxLen = longest;
        size.anyEnds = true;
        size.minLc = c.minLc;
        size.minSize = smallest;
        size.maxSize = c.maxSeparate;

        Step last = new Step();
        int lc = Math.max(0, c.minLc - 10);
        int tolerated = Math.max(PrimerDesign.MASK_TOLERANCE, LAST_MASKED);
        last.relaxed = "Tm " + (c.minTm - 5) + "-" + (c.maxTm + 5) + " °C, primers " + shortest + "-" + longest
                + " nt, any 3' and 5' end, amplicons " + smallest + "-" + c.maxSeparate + " bp, linguistic complexity from "
                + lc + "%, up to " + tolerated + " masked bases, " + LAST_MARGIN + " bp from the primers to the region";
        last.minTm = c.minTm - 5;
        last.maxTm = c.maxTm + 5;
        last.maxTmDifference = 4;
        last.minLen = shortest;
        last.maxLen = longest;
        last.anyEnds = true;
        last.minLc = lc;
        last.maxMasked = tolerated;
        last.minSize = smallest;
        last.maxSize = c.maxSeparate;
        last.margin = LAST_MARGIN;

        return List.of(panel, tm, ends, size, last);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Primers
    // ═══════════════════════════════════════════════════════════════════════
    /** The primer of a step with its 3' end at {@code x}, or null if none passes. */
    private Primer primerAt(boolean forward, int x, Step st) {
        Map<Integer, Primer> known = forward ? st.forward : st.reverse;
        if (!known.containsKey(x)) {
            known.put(x, design(forward, x, st));
        }
        return known.get(x);
    }

    /**
     * Designs the primer whose 3' end is at {@code x}: the end of the binding
     * site (exclusive) for a forward primer, its start for a reverse one. As in
     * the panel, it extends 5' until it reaches the minimum Tm.
     */
    private Primer design(boolean forward, int x, Step st) {
        String strand = forward ? sequence : antisense;
        int end3 = forward ? x : seqLen - x; // exclusive, on the strand of the primer
        if (end3 < st.minLen || end3 > seqLen) {
            return null;
        }
        if (!st.anyEnds && req3 > 0 && !allowed3.contains(strand.substring(end3 - req3, end3))) {
            return null;
        }
        String binding = null;
        double tm = 0;
        for (int len = st.minLen; len <= st.maxLen && len <= end3; len++) {
            String b = strand.substring(end3 - len, end3);
            if (!st.anyEnds && !allowed5.contains(b.substring(0, req5))) {
                continue;
            }
            tm = calculator.getTm(b);
            if (tm >= st.minTm) {
                binding = b;
                break;
            }
        }
        if (binding == null || tm > st.maxTm || !plainBases(binding) || PrimerDesign.containsSSR(binding)) {
            return null;
        }
        Primer p = new Primer();
        int len = binding.length();
        p.start = forward ? x - len : x;
        p.end = p.start + len;
        p.masked = maskedBases(p.start, p.end, st, forward);
        if (p.masked < 0) {
            return null;
        }
        p.lc = calculator.getLC();
        if (p.lc < st.minLc) {
            return null;
        }
        p.gc = calculator.getGC();
        p.tm = tm;
        p.binding = binding;
        p.full = (forward ? c.fTail : c.rTail) + binding;
        if (Oligo.quickDimer(p.full) != 0) {
            return null;
        }
        p.name = forward ? name + ":F_" + (p.start + 1) + "-" + p.end : name + ":R_" + p.end + "-" + (p.start + 1);
        return p;
    }

    private static boolean plainBases(String s) {
        for (int i = 0; i < s.length(); i++) {
            if ("acgt".indexOf(s.charAt(i)) < 0) {
                return false;
            }
        }
        return true;
    }

    /** Masked bases in a binding site, or -1 if too many, or any at its 3' end. */
    private int maskedBases(int start, int end, Step st, boolean forward) {
        int n = 0;
        for (int i = start; i < end; i++) {
            if (masked[i]) {
                n++;
                boolean at3 = forward ? i >= end - CLEAN_3 : i < start + CLEAN_3;
                if (at3 || n > st.maxMasked) {
                    return -1;
                }
            }
        }
        return n;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Pairs
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * The amplicon for what is left uncovered from base {@code pos} on. Fewer
     * PCRs come before stricter filters, since a wider Tm window costs a pair in
     * its own tube nothing: one amplicon for as many stretches as fit, at the
     * first step that allows it, and only then a stretch tiled by several. The
     * last step, whose primers may be less specific, comes after those.
     */
    private Amplicon designAt(int pos, List<int[]> left) {
        List<Step> safe = steps.subList(0, steps.size() - 1);
        int n = 0;
        while (n < left.size() && left.get(n)[1] - pos <= c.maxSeparate) {
            n++;
        }
        for (; n > 0; n--) {
            for (Step st : safe) {
                Amplicon a = amplicon(st, pos, left.get(n - 1)[1], left);
                if (a != null) {
                    return a;
                }
            }
        }
        for (Step st : safe) {
            Amplicon a = amplicon(st, pos, pos + 1, left);
            if (a != null) {
                return a;
            }
        }
        return amplicon(steps.get(steps.size() - 1), pos, pos + 1, left);
    }

    /**
     * The pairs of a step that amplify {@code pos} to {@code end}, or null. The
     * ones amplifying the most of what is left uncovered come first; then the
     * fewest masked bases, the closest Tm (to 1 °C), the shortest amplicon, the
     * highest complexity.
     */
    private Amplicon amplicon(Step st, int pos, int end, List<int[]> left) {
        List<Primer> fs = new ArrayList<>();
        for (int x = Math.max(st.minLen, end - st.maxSize); x <= pos - st.margin; x++) {
            Primer p = primerAt(true, x, st);
            if (p != null) {
                fs.add(p);
            }
        }
        List<Primer> rs = new ArrayList<>();
        for (int x = end + st.margin; x < Math.min(seqLen, pos + st.maxSize); x++) {
            Primer p = primerAt(false, x, st);
            if (p != null) {
                rs.add(p);
            }
        }
        List<Choice> pairs = new ArrayList<>();
        for (Primer f : fs) {
            for (Primer r : rs) {
                int size = r.end - f.start;
                if (size < st.minSize || size > st.maxSize || Math.abs(f.tm - r.tm) > st.maxTmDifference) {
                    continue;
                }
                Choice ch = new Choice();
                ch.f = f;
                ch.r = r;
                ch.lo = f.end + st.margin;
                ch.hi = r.start - st.margin;
                ch.gain = CoverageReport.overlap(left, ch.lo, ch.hi);
                pairs.add(ch);
            }
        }
        pairs.sort(Comparator.comparingInt((Choice p) -> -p.gain)
                .thenComparingInt(p -> p.f.masked + p.r.masked)
                .thenComparingLong(p -> Math.round(p.tmDifference()))
                .thenComparingInt(Choice::size)
                .thenComparingInt(p -> -(p.f.lc + p.r.lc)));
        List<Choice> offered = offer(pairs, end == pos + 1);
        if (offered.isEmpty()) {
            return null;
        }
        Amplicon a = new Amplicon();
        a.step = st;
        a.pairs = offered;
        a.lo = offered.get(0).lo;
        a.hi = offered.get(0).hi;
        for (Choice ch : offered) { // what every pair offered amplifies
            a.lo = Math.max(a.lo, ch.lo);
            a.hi = Math.min(a.hi, ch.hi);
        }
        return a;
    }

    /**
     * True if the 3' end of the primer occurs only at its binding site, on
     * either strand of the target sequence. The masks work on longer k-mers, so
     * a primer can pass them with a 3' end that also fits elsewhere, and in a
     * tube of its own that is where it would misprime. Counted late, for the
     * pairs about to be offered, because it scans the whole sequence.
     */
    private boolean unique3(Primer p) {
        if (p.hits3 < 0) {
            String end = p.binding.substring(Math.max(0, p.binding.length() - UNIQUE_3));
            p.hits3 = occurrences(sequence, end) + occurrences(antisense, end);
        }
        return p.hits3 == 1;
    }

    /** Occurrences of {@code word}, counted up to two: more makes no difference here. */
    private static int occurrences(String text, String word) {
        int first = text.indexOf(word);
        if (first < 0) {
            return 0;
        }
        return (text.indexOf(word, first + 1) < 0) ? 1 : 2;
    }

    /**
     * The best pair free of a primer-dimer, and alternatives with other primers.
     * When tiling, an alternative has to reach as far as the best pair, so that
     * the next amplicon of the tiling fits either.
     */
    private List<Choice> offer(List<Choice> pairs, boolean tiling) {
        List<Choice> out = new ArrayList<>();
        for (Choice ch : pairs) {
            if (out.size() == PAIRS || (tiling && !out.isEmpty() && ch.gain < out.get(0).gain)) {
                break;
            }
            boolean distinct = true;
            for (Choice o : out) {
                distinct &= Math.abs(ch.f.start - o.f.start) >= MIN_SHIFT && Math.abs(ch.r.start - o.r.start) >= MIN_SHIFT;
            }
            if (distinct && unique3(ch.f) && unique3(ch.r) && Oligo.quickDimer(ch.f.full, ch.r.full) == 0) {
                out.add(ch);
            }
        }
        return out;
    }

    private static List<int[]> subtract(List<int[]> stretches, int lo, int hi) {
        List<int[]> out = new ArrayList<>();
        for (int[] v : stretches) {
            if (v[1] <= lo || v[0] >= hi) {
                out.add(v);
                continue;
            }
            if (v[0] < lo) {
                out.add(new int[]{v[0], lo});
            }
            if (v[1] > hi) {
                out.add(new int[]{hi, v[1]});
            }
        }
        return out;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Report
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * Designs the separate PCRs and writes them up; empty when the panel covers
     * every target in full.
     *
     * @param regions the targets not covered in full
     */
    public String report(List<PanelPlot.Uncovered> regions) {
        if (regions.isEmpty()) {
            return "";
        }
        List<int[]> all = new ArrayList<>();
        for (PanelPlot.Uncovered u : regions) {
            all.addAll(u.gaps);
        }
        List<int[]> gaps = CoverageReport.merge(all);

        List<Amplicon> amplicons = new ArrayList<>();
        List<int[]> failed = new ArrayList<>();
        List<int[]> left = gaps;
        while (!left.isEmpty()) {
            int pos = left.get(0)[0];
            Amplicon a = designAt(pos, left);
            if (a != null) {
                amplicons.add(a);
                left = subtract(left, a.lo, a.hi);
            } else { // no pair reaches this base: try further along
                int next = Math.min(left.get(0)[1], pos + SKIP);
                failed.add(new int[]{pos, next});
                left = subtract(left, pos, next);
            }
        }
        failed = CoverageReport.merge(failed);

        StringBuilder s = new StringBuilder(4096);
        s.append("Primer pairs for separate PCR: ").append(name).append("\n\n");
        s.append("The panel does not cover the regions below in full. These pairs were designed without regard to the\n");
        s.append("panel: amplify each amplicon in a tube of its own, at the annealing temperature of its pair, and do not\n");
        s.append("add it to a panel tube. Each amplicon gives its best pair first, then alternatives.\n");
        s.append(genomeChecked
                ? "Repeats in the genome reference were avoided, but only an alignment shows a primer to be specific.\n"
                : "No genome reference was given, so these primers should be checked for specificity before use.\n");

        s.append("\nTarget\tUncovered\tPanel coverage\tWhy the panel has no pair\tAmplicons\n");
        for (PanelPlot.Uncovered u : regions) {
            s.append(u.name).append('\t').append(stretches(u.gaps)).append('\t')
                    .append(String.format(Locale.ROOT, "%.1f%%", 100 * u.coverage)).append('\t')
                    .append(u.why != null ? u.why : "its amplicons do not read this part").append('\t')
                    .append(ampliconsOf(u, amplicons, failed)).append('\n');
        }

        for (int k = 0; k < amplicons.size(); k++) {
            appendAmplicon(s, k + 1, amplicons.get(k), gaps);
        }
        if (!failed.isEmpty()) {
            s.append("\nNo pair was found for ").append(stretches(failed)).append(", even with every filter relaxed.\n");
            s.append("See the panel map for what lies around them, and *_uncovered.fasta for their sequence.\n");
        }
        return s.toString();
    }

    /** The amplicons holding a part of the target, and a note if some of it has no pair. */
    private static String ampliconsOf(PanelPlot.Uncovered u, List<Amplicon> amplicons, List<int[]> failed) {
        List<int[]> merged = CoverageReport.merge(u.gaps);
        List<String> numbers = new ArrayList<>();
        for (int k = 0; k < amplicons.size(); k++) {
            if (CoverageReport.overlap(merged, amplicons.get(k).lo, amplicons.get(k).hi) > 0) {
                numbers.add(String.valueOf(k + 1));
            }
        }
        boolean missing = false;
        for (int[] v : failed) {
            missing |= CoverageReport.overlap(merged, v[0], v[1]) > 0;
        }
        String list = String.join(", ", numbers);
        if (!missing) {
            return list;
        }
        return list.isEmpty() ? "none" : list + ", in part";
    }

    private void appendAmplicon(StringBuilder s, int number, Amplicon a, List<int[]> gaps) {
        List<int[]> holds = new ArrayList<>();
        for (int[] g : gaps) {
            int lo = Math.max(g[0], a.lo);
            int hi = Math.min(g[1], a.hi);
            if (hi > lo) {
                holds.add(new int[]{lo, hi});
            }
        }
        s.append("\nAmplicon ").append(number).append(": covers ").append(stretches(holds)).append('\n');
        s.append("Filters relaxed: ").append(a.step.relaxed).append('\n');
        boolean tails = !c.fTail.isEmpty() || !c.rTail.isEmpty();
        for (int k = 0; k < a.pairs.size(); k++) {
            Choice p = a.pairs.get(k);
            s.append('\n').append(k == 0 ? "Best pair: " : "Alternative " + k + ": ")
                    .append(p.f.start + 1).append('-').append(p.r.end).append(", ").append(p.size()).append(" bp, Tm difference ")
                    .append(String.format(Locale.ROOT, "%.1f", p.tmDifference())).append(" °C\n");
            s.append("Name\tSequence").append(tails ? "\tBinding_sequence" : "").append("\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)\tMasked_bases\n");
            appendPrimer(s, p.f, tails);
            appendPrimer(s, p.r, tails);
            if (p.f.masked + p.r.masked > 0) {
                s.append("Holds repeat-masked bases: check its specificity before use.\n");
            }
            s.append("Amplicon, the uncovered bases in upper case:\n");
            appendAmpliconSequence(s, p, gaps);
        }
    }

    private static void appendPrimer(StringBuilder s, Primer p, boolean tails) {
        s.append(p.name).append('\t').append(p.full);
        if (tails) {
            s.append('\t').append(p.binding);
        }
        s.append('\t').append(p.binding.length()).append('\t').append(String.format(Locale.ROOT, "%.1f", p.tm)).append('\t')
                .append(String.format(Locale.ROOT, "%.1f", p.gc)).append('\t').append(p.lc).append('\t').append(p.masked).append('\n');
    }

    private void appendAmpliconSequence(StringBuilder s, Choice p, List<int[]> gaps) {
        char[] frag = sequence.substring(p.f.start, p.r.end).toCharArray();
        for (int[] g : gaps) {
            for (int i = Math.max(g[0], p.f.start); i < Math.min(g[1], p.r.end); i++) {
                frag[i - p.f.start] = Character.toUpperCase(frag[i - p.f.start]);
            }
        }
        for (int i = 0; i < frag.length; i += LINE) {
            s.append(frag, i, Math.min(LINE, frag.length - i)).append('\n');
        }
    }

    /** 1-based, inclusive coordinates: "133706-133750, 133836-133898". */
    private static String stretches(List<int[]> list) {
        List<String> out = new ArrayList<>();
        for (int[] v : list) {
            out.add((v[0] + 1) + "-" + v[1]);
        }
        return String.join(", ", out);
    }
}
