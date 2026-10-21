import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Panel map (SVG) and BED tracks of one target sequence.
 *
 * The map draws the target regions, coloured by how much of them the panel
 * covers, and under them the amplicons of each tube with their primers. Introns
 * are collapsed: every stretch of sequence that carries a target or an amplicon
 * is drawn to scale, and the stretches between them become a gap labelled with
 * their length, so that the 100 bp exons of a 100 kb gene stay visible.
 *
 * The BED files name the sequence by its record ID, as the primer names do:
 * the primers in the primalscheme layout read by primer-trimming tools, the
 * amplicons coloured by tube, and the targets coloured by coverage.
 */
public class PanelPlot {

    // no yellow or red among the tubes: those mean partly and not covered
    private static final String[] TUBE_COLOURS = {"#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#008300", "#085041", "#854f0b"};
    private static final String COVERED = "#888780";
    private static final String PARTIAL = "#fab219";
    private static final String UNCOVERED = "#d03b3b";
    private static final String INK = "#2c2c2a";
    private static final String MUTED = "#5f5e5a";

    private static final int WIDTH = 1200;
    private static final int LEFT = 90;              // column of track names
    private static final int ROW_W = WIDTH - LEFT - 20;
    private static final double MAX_SCALE = 0.3;     // px per bp
    private static final int GAP = 46;               // px between the stretches of a row
    private static final int TARGET_PAD = 60;        // bp of sequence shown around a target
    private static final int AMPLICON_PAD = 15;
    private static final int JOIN = 100;             // stretches closer than this (bp) are drawn as one
    private static final int HEADER_H = 88;
    private static final int LABELS_H = 30;
    private static final int TARGET_H = 22;
    private static final int LANE_H = 14;
    private static final int NOTE_H = 14;
    private static final int ROW_GAP = 20;

    /** A target region, 0-based half-open. */
    private static final class Target {

        int start;
        int end;
        int number;        // as in the coverage report
        String label;      // e.g. "COL4A5 exon 4"
        String mark;       // drawn above it: exon number, P for a promoter, region number for FASTA
        String gene;
        boolean promoter;
        char strand = '.';
        int widened;       // search flank of its design region, 0 = default
        double coverage;
    }

    /** An amplicon in its tube (1-based). */
    private static final class Amp {

        final Pair pcr;
        final int tube;

        Amp(Pair pcr, int tube) {
            this.pcr = pcr;
            this.tube = tube;
        }

        int start() {
            return pcr.fx5;
        }

        int end() {
            return pcr.rx5 + 1;
        }
    }

    private final String name;
    private final int seqLen;
    private final boolean multiplex;
    private final int ntubes;
    private final List<Target> targets = new ArrayList<>();
    private final List<Amp> amps = new ArrayList<>();
    private final RegionNames names;

    /**
     * @param name      record ID, as in the primer names
     * @param seqLen    length of the target sequence
     * @param bounds    target boundaries before merging, [s1,e1,s2,e2,...]
     * @param merged    the same after merging: a negative start marks a target
     *                  merged into region -start; may be null
     * @param info      annotation of each target, or null for FASTA input
     * @param widened   per design region, the widened search flank (bp), 0 for
     *                  the default; may be null
     * @param tubes     amplicons of this sequence per tube, index 0 = Panel1
     * @param multiplex false for a single-plex panel, drawn as one pool
     */
    public PanelPlot(String name, int seqLen, int[] bounds, int[] merged, List<ExonPromoterExtractor.Region> info,
            int[] widened, List<List<Pair>> tubes, boolean multiplex) {
        this.name = name;
        this.seqLen = seqLen;
        this.multiplex = multiplex;
        this.names = new RegionNames(merged, info);

        for (int k = 0; k < tubes.size(); k++) {
            for (Pair pcr : tubes.get(k)) {
                amps.add(new Amp(pcr, k + 1));
            }
        }
        amps.sort(Comparator.comparingInt(Amp::start).thenComparingInt(Amp::end).thenComparingInt(a -> a.tube));
        int last = multiplex ? 2 : 1;
        for (Amp a : amps) {
            last = Math.max(last, a.tube);
        }
        ntubes = last;

        List<int[]> iv = new ArrayList<>();
        for (Amp a : amps) {
            iv.add(new int[]{a.start(), a.end()});
        }
        List<int[]> covered = CoverageReport.merge(iv);

        int n = 0;
        String previousGene = null;
        for (int k = 0; k + 1 < bounds.length; k += 2) {
            int[] t0 = CoverageReport.targetInterval(bounds[k], bounds[k + 1]);
            if (t0 == null) {
                continue;
            }
            n++;
            Target t = new Target();
            t.start = Math.max(0, t0[0]);
            t.end = Math.min(seqLen, t0[1]);
            t.number = n;
            ExonPromoterExtractor.Region r = (info != null && k / 2 < info.size()) ? info.get(k / 2) : null;
            if (r != null) {
                t.gene = r.gene;
                t.promoter = "PROMOTER".equals(r.type);
                t.strand = r.strand;
                t.label = RegionNames.label(r);
                t.mark = t.promoter ? "P" : r.idx;
                if (!r.gene.equals(previousGene)) {
                    t.mark = r.gene + " " + t.mark; // name the gene where its regions begin
                }
                previousGene = r.gene;
            } else {
                t.label = "region " + n;
                t.mark = String.valueOf(n);
            }
            int design = (merged != null && k < merged.length && merged[k] < 0) ? -merged[k] - 1 : k / 2;
            t.widened = (widened != null && design < widened.length) ? widened[design] : 0;
            t.coverage = (t.end > t.start) ? CoverageReport.overlap(covered, t.start, t.end) / (double) (t.end - t.start) : 1;
            targets.add(t);
        }
    }

    private String tubeName(int tube) {
        return multiplex ? "Panel" + tube : "Panel";
    }

    private static String tubeColour(int tube) {
        return TUBE_COLOURS[(tube - 1) % TUBE_COLOURS.length];
    }

    private static String status(Target t) {
        return (t.coverage >= 1) ? COVERED : (t.coverage > 0 ? PARTIAL : UNCOVERED);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // SVG
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * The panel map as a standalone SVG document; hovering over a region or an
     * amplicon shows its details.
     */
    public String svg() {
        // Stretches of sequence drawn to scale: targets and amplicons with some context
        List<int[]> raw = new ArrayList<>();
        for (Target t : targets) {
            raw.add(new int[]{Math.max(0, t.start - TARGET_PAD), Math.min(seqLen, t.end + TARGET_PAD)});
        }
        for (Amp a : amps) {
            raw.add(new int[]{Math.max(0, a.start() - AMPLICON_PAD), Math.min(seqLen, a.end() + AMPLICON_PAD)});
        }
        raw.sort(Comparator.comparingInt(v -> v[0]));
        List<int[]> st = new ArrayList<>(); // {lo, hi}
        for (int[] v : raw) {
            int[] prev = st.isEmpty() ? null : st.get(st.size() - 1);
            if (prev != null && v[0] <= prev[1] + JOIN) {
                prev[1] = Math.max(prev[1], v[1]);
            } else {
                st.add(new int[]{v[0], v[1]});
            }
        }
        int longest = 1;
        for (int[] s : st) {
            longest = Math.max(longest, s[1] - s[0]);
        }
        double scale = Math.min(MAX_SCALE, ROW_W / (double) longest);

        // Rows: stretches are laid out left to right and wrapped
        int ns = st.size();
        double[] sx = new double[ns];
        int[] srow = new int[ns];
        double x = 0;
        int rows = 1;
        for (int i = 0; i < ns; i++) {
            double w = (st.get(i)[1] - st.get(i)[0]) * scale;
            if (x > 0 && x + w > ROW_W) {
                rows++;
                x = 0;
            }
            sx[i] = x;
            srow[i] = rows - 1;
            x += w + GAP;
        }

        // Lanes: amplicons of one tube that overlap on screen go to separate lanes
        List<List<List<double[]>>> laneEnds = new ArrayList<>(); // [row][tube] -> lane ends
        for (int r = 0; r < rows; r++) {
            List<List<double[]>> perTube = new ArrayList<>();
            for (int t = 0; t <= ntubes; t++) {
                perTube.add(new ArrayList<>());
            }
            laneEnds.add(perTube);
        }
        int[] lane = new int[amps.size()];
        for (int i = 0; i < amps.size(); i++) {
            Amp a = amps.get(i);
            int s = stretchOf(st, a.start());
            double x1 = sx[s] + (a.start() - st.get(s)[0]) * scale;
            double x2 = sx[s] + (a.end() - st.get(s)[0]) * scale;
            List<double[]> ends = laneEnds.get(srow[s]).get(a.tube);
            int l = 0;
            while (l < ends.size() && ends.get(l)[0] > x1 - 2) {
                l++;
            }
            if (l == ends.size()) {
                ends.add(new double[]{x2});
            } else {
                ends.get(l)[0] = x2;
            }
            lane[i] = l;
        }

        // Vertical layout of each row
        double[] rowY = new double[rows];
        double[][] tubeY = new double[rows][ntubes + 1];
        double[] noteY = new double[rows];
        double y = HEADER_H;
        for (int r = 0; r < rows; r++) {
            rowY[r] = y;
            double ty = y + LABELS_H + TARGET_H + 6;
            for (int t = 1; t <= ntubes; t++) {
                tubeY[r][t] = ty;
                ty += Math.max(1, laneEnds.get(r).get(t).size()) * LANE_H + 4;
            }
            noteY[r] = ty + 8;
            y = ty + NOTE_H + ROW_GAP;
        }
        int height = (int) Math.ceil(y);

        StringBuilder s = new StringBuilder(1 << 16);
        s.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        s.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"").append(WIDTH).append("\" height=\"").append(height)
                .append("\" viewBox=\"0 0 ").append(WIDTH).append(' ').append(height)
                .append("\" font-family=\"Arial, Helvetica, sans-serif\">\n");
        s.append("<title>").append(esc(name)).append(" panel map</title>\n");
        s.append("<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n");
        header(s);

        for (int r = 0; r < rows; r++) {
            double targetMid = rowY[r] + LABELS_H + TARGET_H / 2.0;
            text(s, 12, targetMid, "Targets", 12, MUTED, "start");
            for (int t = 1; t <= ntubes; t++) {
                text(s, 12, tubeY[r][t] + 5, tubeName(t), 12, MUTED, "start");
            }
        }

        // Backbone of each stretch, and the length of the sequence left out between stretches
        for (int i = 0; i < ns; i++) {
            double mid = rowY[srow[i]] + LABELS_H + TARGET_H / 2.0;
            double x1 = LEFT + sx[i];
            double x2 = x1 + (st.get(i)[1] - st.get(i)[0]) * scale;
            s.append("<line x1=\"").append(f(x1)).append("\" y1=\"").append(f(mid)).append("\" x2=\"").append(f(x2))
                    .append("\" y2=\"").append(f(mid)).append("\" stroke=\"#b4b2a9\" stroke-width=\"1\"/>\n");
            if (i + 1 < ns && srow[i + 1] == srow[i]) {
                text(s, x2 + GAP / 2.0, mid, gapLength(st.get(i + 1)[0] - st.get(i)[1]), 10, "#888780", "middle");
            }
        }

        // Target regions, their marks above and, when not covered in full, a note below
        List<List<double[]>> marks = new ArrayList<>(); // per row: occupied spans on the two mark lines and the note line
        for (int r = 0; r < rows; r++) {
            marks.add(new ArrayList<>(List.of(new double[]{-1e9}, new double[]{-1e9}, new double[]{-1e9})));
        }
        for (Target t : targets) {
            int i = stretchOf(st, t.start);
            int r = srow[i];
            double x1 = LEFT + sx[i] + (t.start - st.get(i)[0]) * scale;
            double w = Math.max(2, (t.end - t.start) * scale);
            double ty = rowY[r] + LABELS_H + (TARGET_H - 14) / 2.0;
            String colour = status(t);
            s.append("<rect x=\"").append(f(x1)).append("\" y=\"").append(f(ty)).append("\" width=\"").append(f(w))
                    .append("\" height=\"14\" rx=\"2\" fill=\"").append(colour).append('"');
            if (t.promoter) {
                s.append(" fill-opacity=\"0.3\" stroke=\"").append(colour).append("\" stroke-width=\"1.2\"");
            }
            s.append("><title>").append(esc(targetTip(t))).append("</title></rect>\n");

            double cx = x1 + w / 2;
            double mw = t.mark.length() * 6.0 + 4;
            List<double[]> m = marks.get(r);
            if (cx - mw / 2 > m.get(0)[0]) {
                text(s, cx, rowY[r] + 22, t.mark, 11, INK, "middle");
                m.get(0)[0] = cx + mw / 2;
            } else if (cx - mw / 2 > m.get(1)[0]) {
                text(s, cx, rowY[r] + 9, t.mark, 11, INK, "middle");
                m.get(1)[0] = cx + mw / 2;
            }
            if (t.coverage < 1) {
                String note = (t.coverage > 0) ? pct(t.coverage) : "not covered";
                double nw = note.length() * 5.6 + 4;
                if (cx - nw / 2 > m.get(2)[0]) {
                    text(s, cx, noteY[r], note, 10, (t.coverage > 0) ? "#ba7517" : UNCOVERED, "middle");
                    m.get(2)[0] = cx + nw / 2;
                }
            }
        }

        // Amplicons: the amplicon pale, its primers in full colour
        for (int i = 0; i < amps.size(); i++) {
            Amp a = amps.get(i);
            Pair p = a.pcr;
            int si = stretchOf(st, a.start());
            int r = srow[si];
            double x1 = LEFT + sx[si] + (a.start() - st.get(si)[0]) * scale;
            double x2 = LEFT + sx[si] + (a.end() - st.get(si)[0]) * scale;
            double ay = tubeY[r][a.tube] + lane[i] * LANE_H;
            String colour = tubeColour(a.tube);
            s.append("<g><title>").append(esc(ampliconTip(a))).append("</title>");
            s.append("<rect x=\"").append(f(x1)).append("\" y=\"").append(f(ay)).append("\" width=\"").append(f(x2 - x1))
                    .append("\" height=\"10\" rx=\"2\" fill=\"").append(colour).append("\" fill-opacity=\"0.35\"");
            if (p.note != null) {
                s.append(" stroke=\"").append(UNCOVERED).append("\" stroke-dasharray=\"3 2\"");
            }
            s.append("/>");
            s.append("<rect x=\"").append(f(x1)).append("\" y=\"").append(f(ay)).append("\" width=\"").append(f(Math.max(1, p.fln * scale)))
                    .append("\" height=\"10\" fill=\"").append(colour).append("\"/>");
            double rw = Math.max(1, p.rln * scale);
            s.append("<rect x=\"").append(f(x2 - rw)).append("\" y=\"").append(f(ay)).append("\" width=\"").append(f(rw))
                    .append("\" height=\"10\" fill=\"").append(colour).append("\"/></g>\n");
        }
        s.append("</svg>\n");
        return s.toString();
    }

    /** Title, summary line and legend. */
    private void header(StringBuilder s) {
        text(s, 12, 26, name + " — panel map", 16, INK, "start");
        int targetBp = 0;
        double coveredBp = 0;
        for (Target t : targets) {
            targetBp += t.end - t.start;
            coveredBp += t.coverage * (t.end - t.start);
        }
        int[] perTube = new int[ntubes + 1];
        for (Amp a : amps) {
            perTube[a.tube]++;
        }
        StringBuilder sum = new StringBuilder();
        sum.append(targets.size()).append(" target region(s), ").append(targetBp).append(" bp; covered ")
                .append(targetBp > 0 ? pct(coveredBp / targetBp) : "n/a").append("; ");
        for (int t = 1; t <= ntubes; t++) {
            sum.append(t > 1 ? ", " : "").append(tubeName(t)).append(' ').append(perTube[t]);
        }
        sum.append(" amplicons. Introns collapsed; hover over a region or amplicon for details.");
        text(s, 12, 46, sum.toString(), 12, MUTED, "start");

        double x = 12;
        double y = 68;
        for (int t = 1; t <= ntubes; t++) {
            x = legend(s, x, y, tubeColour(t), tubeName(t), false);
        }
        x += 14;
        x = legend(s, x, y, COVERED, "covered", false);
        x = legend(s, x, y, PARTIAL, "partly covered", false);
        x = legend(s, x, y, UNCOVERED, "not covered", false);
        boolean promoters = false;
        for (Target t : targets) {
            promoters |= t.promoter;
        }
        if (promoters) {
            x = legend(s, x, y, COVERED, "promoter", true);
        }
        text(s, x + 10, y, "Dark ends of an amplicon: its primers", 11, MUTED, "start");
    }

    private double legend(StringBuilder s, double x, double y, String colour, String label, boolean outlined) {
        s.append("<rect x=\"").append(f(x)).append("\" y=\"").append(f(y - 5)).append("\" width=\"10\" height=\"10\" rx=\"2\" fill=\"")
                .append(colour).append('"');
        if (outlined) {
            s.append(" fill-opacity=\"0.3\" stroke=\"").append(colour).append("\" stroke-width=\"1.2\"");
        }
        s.append("/>\n");
        text(s, x + 14, y, label, 11, MUTED, "start");
        return x + 14 + label.length() * 6.0 + 14;
    }

    private String targetTip(Target t) {
        StringBuilder s = new StringBuilder();
        s.append(t.label).append(" (region ").append(t.number).append(")\n");
        s.append(t.start + 1).append('-').append(t.end).append(", ").append(t.end - t.start).append(" bp");
        if (t.strand == '+' || t.strand == '-') {
            s.append(", ").append(t.strand).append(" strand");
        }
        s.append("\nCovered: ").append(pct(t.coverage));
        if (t.coverage <= 0) {
            s.append(" (no primer pair reaches it)");
        }
        if (t.widened > 0) {
            s.append("\nSearch flanks widened to ").append(t.widened).append(" bp");
        }
        return s.toString();
    }

    private String ampliconTip(Amp a) {
        Pair p = a.pcr;
        StringBuilder s = new StringBuilder();
        s.append(tubeName(a.tube)).append(", ").append(names.design(p.exonid)).append('\n');
        s.append(a.start() + 1).append('-').append(a.end()).append(", ").append(p.pcrsize).append(" bp\n");
        s.append(p.fprimername).append(" (").append(p.fln).append(" nt)\n");
        s.append(p.rprimername).append(" (").append(p.rln).append(" nt)");
        if (p.note != null) {
            s.append("\nFORCED: ").append(p.note);
        }
        return s.toString();
    }

    /** Index of the stretch that holds position {@code b}. */
    private static int stretchOf(List<int[]> st, int b) {
        int lo = 0;
        int hi = st.size() - 1;
        while (lo < hi) {
            int mid = (lo + hi + 1) >>> 1;
            if (st.get(mid)[0] <= b) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        return lo;
    }

    private static String gapLength(int bp) {
        return (bp >= 1000) ? String.format("%.1fk", bp / 1000.0) : bp + " bp";
    }

    private static String pct(double v) {
        return String.format("%.1f%%", 100 * v);
    }

    private static void text(StringBuilder s, double x, double y, String text, int size, String colour, String anchor) {
        s.append("<text x=\"").append(f(x)).append("\" y=\"").append(f(y)).append("\" font-size=\"").append(size)
                .append("\" fill=\"").append(colour).append("\" text-anchor=\"").append(anchor)
                .append("\" dominant-baseline=\"central\">").append(esc(text)).append("</text>\n");
    }

    /** SVG coordinates always use a decimal point, whatever the locale. */
    private static String f(double v) {
        return String.format(Locale.ROOT, "%.1f", v);
    }

    private static String esc(String v) {
        return v.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // BED
    // ═══════════════════════════════════════════════════════════════════════
    /**
     * Panel primers in the primalscheme BED layout: chrom, start, end,
     * name_N_LEFT / name_N_RIGHT, tube (pool), strand, primer sequence without
     * the adapter tail. Amplicons are numbered along the sequence. This is the
     * layout read by primer-trimming tools such as ivar trim and samtools
     * ampliconclip.
     */
    public String primersBed() {
        StringBuilder s = new StringBuilder();
        int n = 0;
        for (Amp a : amps) {
            n++;
            Pair p = a.pcr;
            s.append(name).append('\t').append(p.fx5).append('\t').append(p.fx5 + p.fln).append('\t')
                    .append(name).append('_').append(n).append("_LEFT\t").append(a.tube).append("\t+\t")
                    .append(binding(p.fprimer, p.fln)).append('\n');
            s.append(name).append('\t').append(p.rx5 + 1 - p.rln).append('\t').append(p.rx5 + 1).append('\t')
                    .append(name).append('_').append(n).append("_RIGHT\t").append(a.tube).append("\t-\t")
                    .append(binding(p.rprimer, p.rln)).append('\n');
        }
        return s.toString();
    }

    /**
     * Amplicons for a genome browser, coloured by tube (also the score column);
     * the thick part is the insert between the primers.
     */
    public String ampliconsBed() {
        StringBuilder s = new StringBuilder();
        s.append("track name=\"").append(name).append(" amplicons\" description=\"PCRpanel amplicons by tube (score = tube); thin ends are the primers\" itemRgb=\"On\"\n");
        int n = 0;
        for (Amp a : amps) {
            n++;
            Pair p = a.pcr;
            s.append(name).append('\t').append(a.start()).append('\t').append(a.end()).append('\t')
                    .append(name).append('_').append(n).append('\t').append(a.tube).append("\t.\t")
                    .append(a.start() + p.fln).append('\t').append(a.end() - p.rln).append('\t')
                    .append(rgb(tubeColour(a.tube))).append('\n');
        }
        return s.toString();
    }

    /**
     * Target regions for a genome browser, coloured by coverage; the score
     * column is the covered fraction in thousandths.
     */
    public String targetsBed() {
        StringBuilder s = new StringBuilder();
        s.append("track name=\"").append(name).append(" targets\" description=\"PCRpanel targets: grey covered, amber partly covered, red not covered\" itemRgb=\"On\"\n");
        for (Target t : targets) {
            s.append(name).append('\t').append(t.start).append('\t').append(t.end).append('\t')
                    .append("region").append(t.number).append('_').append(t.label.replaceAll("[^A-Za-z0-9_.-]+", "_")).append('\t')
                    .append(Math.round(t.coverage * 1000)).append('\t').append(t.strand).append('\t')
                    .append(t.start).append('\t').append(t.end).append('\t').append(rgb(status(t))).append('\n');
        }
        return s.toString();
    }

    /** The part of a primer that binds the target: its 3' end, after the tail. */
    private static String binding(String primer, int length) {
        return primer.substring(Math.max(0, primer.length() - length));
    }

    private static String rgb(String hex) {
        return Integer.parseInt(hex.substring(1, 3), 16) + "," + Integer.parseInt(hex.substring(3, 5), 16) + ","
                + Integer.parseInt(hex.substring(5, 7), 16);
    }
}
