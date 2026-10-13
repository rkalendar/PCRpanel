import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Target-coverage report for a designed panel.
 *
 * Amplicon count is not a measure of coverage, because amplicons overlap and
 * because regions excluded by repeat or low-complexity masking are silently
 * skipped. This class measures, base by base, how much of the annotated target
 * each pool actually covers, and how much the pools cover together, so that any
 * gap left by the design is visible in the output rather than having to be
 * reconstructed from the amplicon coordinates afterwards.
 *
 * Coordinates are half-open [start, end) on the target sequence, matching the
 * primer coordinates used elsewhere in the tool.
 */
public class CoverageReport {

    /** One annotated target region (an exon, or a merged group of exons). */
    private static final class Region {

        final int start;
        final int end;
        final int index;

        Region(int index, int start, int end) {
            this.index = index;
            this.start = start;
            this.end = end;
        }

        int length() {
            return end - start;
        }
    }

    private final List<Region> regions = new ArrayList<>();
    private final Map<String, List<int[]>> pools = new LinkedHashMap<>();

    /**
     * @param exonBoundaries pairs of start/end positions, as passed to the
     *                       designer before exon merging. Negative entries mark
     *                       regions folded into a neighbour and are skipped.
     */
    public CoverageReport(int[] exonBoundaries) {
        int n = 0;
        for (int i = 0; i + 1 < exonBoundaries.length; i += 2) {
            int[] t = targetInterval(exonBoundaries[i], exonBoundaries[i + 1]);
            if (t == null) {
                continue;
            }
            n++;
            regions.add(new Region(n, t[0], t[1]));
        }
    }

    /**
     * One target region in the amplicon frame, or null for a negative (skipped)
     * entry. Exon boundaries arrive 1-based inclusive [s, e], while amplicon
     * coordinates are 0-based half-open (fx5 is a 0-based offset, the printed
     * primer label adds 1).
     */
    static int[] targetInterval(int s, int e) {
        if (s < 0 || e < 0) {
            return null;
        }
        if (e < s) {
            int t = s;
            s = e;
            e = t;
        }
        return new int[]{s - 1, e};
    }

    /** Register one pool and the amplicons it carries. */
    public void addPool(String label, List<int[]> amplicons) {
        pools.put(label, amplicons);
    }

    public boolean isEmpty() {
        return regions.isEmpty() || pools.isEmpty();
    }

    static List<int[]> merge(List<int[]> intervals) {
        List<int[]> sorted = new ArrayList<>(intervals);
        sorted.sort((a, b) -> Integer.compare(a[0], b[0]));
        List<int[]> out = new ArrayList<>();
        for (int[] iv : sorted) {
            if (!out.isEmpty() && iv[0] <= out.get(out.size() - 1)[1]) {
                int[] last = out.get(out.size() - 1);
                last[1] = Math.max(last[1], iv[1]);
            } else {
                out.add(new int[]{iv[0], iv[1]});
            }
        }
        return out;
    }

    static int overlap(List<int[]> merged, int start, int end) {
        int total = 0;
        for (int[] iv : merged) {
            int lo = Math.max(start, iv[0]);
            int hi = Math.min(end, iv[1]);
            if (hi > lo) {
                total += hi - lo;
            }
        }
        return total;
    }

    private int targetBases() {
        int t = 0;
        for (Region r : regions) {
            t += r.length();
        }
        return t;
    }

    private static String pct(int part, int whole) {
        return whole == 0 ? "n/a" : String.format("%.1f%%", 100.0 * part / whole);
    }

    /** Human-readable report, appended to the panel output file. */
    public String render() {
        if (isEmpty()) {
            return "";
        }
        StringBuilder s = new StringBuilder(4096);
        int targetBp = targetBases();

        s.append("\n\nTarget coverage report\n");
        s.append("Annotated target: ").append(regions.size()).append(" region(s), ")
                .append(targetBp).append(" bp\n\n");
        s.append("Pool\tAmplicons\tTarget bp covered\tCoverage\tRegions full\tPartial\tUncovered\n");

        List<int[]> all = new ArrayList<>();
        for (Map.Entry<String, List<int[]>> e : pools.entrySet()) {
            all.addAll(e.getValue());
            s.append(summarise(e.getKey(), e.getValue(), targetBp));
        }
        if (pools.size() > 1) {
            s.append(summarise("All pools combined", all, targetBp));
        }

        List<int[]> merged = merge(all);
        StringBuilder gaps = new StringBuilder();
        int nGaps = 0;
        int gapBp = 0;
        for (Region r : regions) {
            int c = overlap(merged, r.start, r.end);
            if (c < r.length()) {
                nGaps++;
                gapBp += r.length() - c;
                if (gaps.length() > 0) {
                    gaps.append(", ");
                }
                gaps.append(r.index).append(" (").append(pct(c, r.length())).append(")");
            }
        }
        if (nGaps == 0) {
            s.append("\nAll annotated target regions are covered in full.\n");
        } else {
            s.append("\nRegions not covered in full: ").append(nGaps).append(" of ")
                    .append(regions.size()).append(", leaving ").append(gapBp)
                    .append(" bp of the target without an amplicon.\n");
            s.append("Region (coverage): ").append(gaps).append("\n");
            s.append("\nA region is left uncovered when no primer pair satisfying the current length, Tm,\n");
            s.append("complexity and masking constraints could be placed around it. To close such gaps,\n");
            s.append("widen minPCR/maxPCR, minLen/maxLen and minTm/maxTm, or relax the masking settings\n");
            s.append("(minLC, maskedBases, repeatKmer, repeatMinCount), and re-run; primers recovered this\n");
            s.append("way are less stringently filtered and should be checked for specificity before use.\n");
        }
        return s.toString();
    }

    private String summarise(String label, List<int[]> amplicons, int targetBp) {
        List<int[]> merged = merge(amplicons);
        int covered = 0;
        int full = 0;
        int partial = 0;
        int none = 0;
        for (Region r : regions) {
            int c = overlap(merged, r.start, r.end);
            covered += c;
            if (c == 0) {
                none++;
            } else if (c >= r.length()) {
                full++;
            } else {
                partial++;
            }
        }
        return label + "\t" + amplicons.size() + "\t" + covered + "\t"
                + pct(covered, targetBp) + "\t" + full + "\t" + partial + "\t" + none + "\n";
    }

    /** Convenience: amplicon interval from a designed pair. */
    public static int[] interval(int start, int end) {
        return new int[]{Math.min(start, end), Math.max(start, end)};
    }
}
