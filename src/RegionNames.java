import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Names of the target regions of one sequence, as the reports and the panel
 * map show them: gene and exon number from the annotation, then the region
 * number, e.g. "COL4A5 exon 5 (region 11)".
 *
 * Region numbers count every target of the sequence in order, promoters and
 * the exons of neighbouring genes included, so they differ from the exon
 * numbers of the annotation; the name carries both. Without an annotation
 * (FASTA input) a region is named by its number alone.
 */
public class RegionNames {

    private final int[] merged;
    private final List<ExonPromoterExtractor.Region> info;

    /**
     * @param merged target boundaries after merging, [s1,e1,s2,e2,...]: a
     *               negative start marks a target merged into region -start;
     *               may be null
     * @param info   annotation of each boundary pair, or null for FASTA input
     */
    public RegionNames(int[] merged, List<ExonPromoterExtractor.Region> info) {
        this.merged = merged;
        this.info = info;
    }

    /** One annotated target, e.g. "COL4A5 exon 4" or "COL4A5 promoter". */
    static String label(ExonPromoterExtractor.Region r) {
        return r.gene + (isPromoter(r) ? " promoter" : " exon " + r.idx);
    }

    private static boolean isPromoter(ExonPromoterExtractor.Region r) {
        return "PROMOTER".equals(r.type);
    }

    private ExonPromoterExtractor.Region annotation(int k) {
        return (info != null && k >= 0 && k < info.size()) ? info.get(k) : null;
    }

    /** The design region (0-based) that target {@code k} was merged into, or {@code k} itself. */
    public int designOf(int k) {
        return (merged != null && 2 * k < merged.length && merged[2 * k] < 0) ? -merged[2 * k] - 1 : k;
    }

    /** Target {@code k} (0-based) on its own, e.g. "COL4A5 exon 6 (region 12)". */
    public String target(int k) {
        ExonPromoterExtractor.Region r = annotation(k);
        return (r != null) ? label(r) + " (region " + (k + 1) + ")" : "region " + (k + 1);
    }

    /**
     * Design region {@code h} (0-based, as {@link Pair#exonid}) with the targets
     * merged into it, e.g. "COL4A5 exons 5-6 (region 11)".
     */
    public String design(int h) {
        List<ExonPromoterExtractor.Region> group = new ArrayList<>();
        for (int k = 0; info != null && k < info.size(); k++) {
            boolean member = (k == h) || (merged != null && 2 * k < merged.length && merged[2 * k] == -(h + 1));
            if (member) {
                group.add(info.get(k));
            }
        }
        if (group.isEmpty() || annotation(h) == null) {
            return "region " + (h + 1);
        }
        return summary(group) + " (region " + (h + 1) + ")";
    }

    /**
     * Per gene, in order of appearance, its promoter and its exons with runs of
     * numbers joined: "COL4A6 promoter, exons 1-3; COL4A5 exon 1".
     */
    private static String summary(List<ExonPromoterExtractor.Region> group) {
        Map<String, Boolean> promoter = new LinkedHashMap<>();
        Map<String, TreeSet<Integer>> numbers = new LinkedHashMap<>();
        Map<String, Set<String>> others = new LinkedHashMap<>(); // exon numbers that are not plain integers
        for (ExonPromoterExtractor.Region r : group) {
            promoter.merge(r.gene, isPromoter(r), Boolean::logicalOr);
            numbers.putIfAbsent(r.gene, new TreeSet<>());
            others.putIfAbsent(r.gene, new LinkedHashSet<>());
            if (!isPromoter(r)) {
                try {
                    numbers.get(r.gene).add(Integer.parseInt(r.idx.trim()));
                } catch (NumberFormatException e) {
                    others.get(r.gene).add(r.idx);
                }
            }
        }
        StringBuilder s = new StringBuilder();
        for (String gene : promoter.keySet()) {
            List<String> parts = new ArrayList<>();
            if (promoter.get(gene)) {
                parts.add("promoter");
            }
            List<String> exons = runs(numbers.get(gene));
            exons.addAll(others.get(gene));
            if (!exons.isEmpty()) {
                int n = numbers.get(gene).size() + others.get(gene).size();
                parts.add((n == 1 ? "exon " : "exons ") + String.join(", ", exons));
            }
            s.append(s.length() > 0 ? "; " : "").append(gene).append(' ').append(String.join(", ", parts));
        }
        return s.toString();
    }

    /** 1, 2, 3, 5 -> "1-3", "5". */
    private static List<String> runs(TreeSet<Integer> numbers) {
        List<String> out = new ArrayList<>();
        Integer first = null;
        Integer last = null;
        for (int n : numbers) {
            if (last != null && n == last + 1) {
                last = n;
                continue;
            }
            if (first != null) {
                out.add(first.equals(last) ? String.valueOf(first) : first + "-" + last);
            }
            first = n;
            last = n;
        }
        if (first != null) {
            out.add(first.equals(last) ? String.valueOf(first) : first + "-" + last);
        }
        return out;
    }
}
