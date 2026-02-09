import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

public class NcRNAExonPromoterExtractor {

    private final Path gbPath;
    private final int promoterLen;
    private final List<Region> regions;   // <-- what you need to "return"
    private final String sequence;        // DNA from ORIGIN    
    private final String accession;   // например: NG_034200
    private final String version;     // например: NG_034200.2
    private final String recordId;    // version если есть, иначе accession

    public NcRNAExonPromoterExtractor(String arg) throws Exception {
        this(arg, 1000);
    }

    public NcRNAExonPromoterExtractor(String arg, int promoterLen) throws Exception {
        if (arg == null || arg.isBlank()) {
            throw new IllegalArgumentException("GenBank path is empty");
        }
        this.gbPath = Paths.get(arg);
        Header h = readHeaderFromGenBank(this.gbPath);
        this.accession = h.accession;
        this.version = h.version;
        this.recordId = (h.version != null && !h.version.isBlank()) ? h.version : h.accession;

        this.promoterLen = promoterLen;

        // Save the calculated regions in the field
        this.regions = Collections.unmodifiableList(extract(this.gbPath, this.promoterLen));
        this.sequence = readSequenceFromGenBank(this.gbPath);
    }

    static final class Header {
        final String accession;
        final String version;

        Header(String accession, String version) {
            this.accession = accession == null ? "" : accession;
            this.version = version == null ? "" : version;
        }
    }

    /**
     *
     * @param genBankFile
     * @return
     * @throws IOException
     */
    public static Header readHeaderFromGenBank(Path genBankFile) throws IOException {
        String acc = "";
        String ver = "";

        try (BufferedReader br = Files.newBufferedReader(genBankFile, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                String t = line.trim();

                // It is usually sufficient to read up to FEATURES/ORIGIN.
                if (t.startsWith("FEATURES") || t.startsWith("ORIGIN")) {
                    break;
                }

                if (t.startsWith("ACCESSION")) {
                    // ACCESSION   NG_034200
                    String rest = t.substring("ACCESSION".length()).trim();
                    if (!rest.isEmpty()) {
                        acc = rest.split("\\s+")[0]; // берём первый (primary) accession
                    }
                } else if (t.startsWith("VERSION")) {
                    // VERSION     NG_034200.2  GI:...
                    String rest = t.substring("VERSION".length()).trim();
                    if (!rest.isEmpty()) {
                        ver = rest.split("\\s+")[0]; // первый токен = accession.version
                    }
                }

                if (!acc.isEmpty() && !ver.isEmpty()) {
                    break;
                }
            }
        }

        return new Header(acc, ver);
    }

    public List<Region> getRegions() {
        return regions;
    }

    // optional: move printing to a method rather than the constructor
    public void printTSV() {
        System.out.println("transcript_id\tgene\ttype\tidx\tstart\tend\tstrand");
        for (Region r : regions) {
            System.out.printf("%s\t%s\t%s\t%s\t%d\t%d\t%c%n",
                    nvl(r.transcriptId), nvl(r.gene), r.type, nvl(r.idx),
                    r.start, r.end, r.strand);
        }
    }

    public String getSequence() {
        return sequence;
    }

    public String getAccession() {
        return accession;
    }

    public String getVersion() {
        return version;
    }

    public String getRecordId() {
        return recordId;
    } 

    public int getSequenceLength() {
        return sequence.length();
    }

    public static String readSequenceFromGenBank(Path genBankFile) throws IOException {
        Objects.requireNonNull(genBankFile, "genBankFile");

        StringBuilder sb = new StringBuilder(1024 * 16);
        boolean inOrigin = false;

        try (BufferedReader br = Files.newBufferedReader(genBankFile, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {

                if (!inOrigin) { // robust: ORIGIN / origin / "  ORIGIN"                   
                    String t = line.trim();
                    if (t.regionMatches(true, 0, "ORIGIN", 0, 6)) {
                        inOrigin = true;
                    }
                    continue;
                }

                // end of sequence section
                if (line.startsWith("//")) {
                    break;
                }

                // GenBank ORIGIN lines: digits + spaces + sequence blocks
                for (int i = 0; i < line.length(); i++) {
                    char c = line.charAt(i);
                    if (c >= 'A' && c <= 'Z') {
                        char u = (char) (c - 32); // toLowerCase for ASCII
                        if (u == 'U') {
                            u = 't';
                        }
                        if (u == 'I') {
                            u = 'g';
                        }
                        sb.append(u);
                    } else if (c >= 'a' && c <= 'z') {
                        char u = c;
                        if (u == 'U') {
                            u = 't';
                        }
                        if (u == 'i') {
                            u = 'g';
                        }
                        sb.append(u);
                    }
                }
            }
        }
        return dna.DNA(sb.toString());
    }

    /**
     * Extract regions from a GenBank file.
     *
     * @param genBankFile
     * @param promoterLen
     * @return
     * @throws java.io.IOException
     */
    public static List<Region> extract(Path genBankFile, int promoterLen) throws IOException {
        Objects.requireNonNull(genBankFile, "genBankFile");

        // mutable box for seqLen inference (LOCUS or source feature)
        final int[] seqLenBox = new int[]{-1};

        boolean inFeatures = false;
        Feature cur = null;
        List<Region> out = new ArrayList<>();

        try (BufferedReader br = Files.newBufferedReader(genBankFile, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {

                // 1) try parse length from LOCUS
                if (seqLenBox[0] < 0 && line.startsWith("LOCUS")) {
                    Integer len = parseLocusLength(line);
                    if (len != null) {
                        seqLenBox[0] = len;
                    }
                }

                if (line.startsWith("FEATURES")) {
                    inFeatures = true;
                    continue;
                }
                if (!inFeatures) {
                    continue;
                }

                if (line.startsWith("ORIGIN")) {
                    if (cur != null) {
                        processFeature(cur, seqLenBox, promoterLen, out);
                    }
                    break;
                }

                if (line.length() >= 21 && line.startsWith("     ")) {
                    String key = safeSubstring(line, 5, 21).trim();
                    String rest = safeSubstring(line, 21, line.length()).trim();

                    if (!key.isEmpty()) {
                        // new feature begins
                        if (cur != null) {
                            processFeature(cur, seqLenBox, promoterLen, out);
                        }
                        cur = new Feature(key, rest);
                    } else if (cur != null) {
                        // continuation line (location or qualifiers)
                        cur.consumeContinuation(rest);
                    }
                } else if (cur != null && line.startsWith("                     ")) {
                    String rest = safeSubstring(line, 21, line.length()).trim();
                    cur.consumeContinuation(rest);
                }
            }
        }

        return mergeAndDedupByGene(out, false); // false: не сливать "впритык"
    }

    static List<Region> mergeAndDedupByGene(List<Region> regions, boolean mergeTouching) {
        // Group by (gene, type, strand) — this removes overlaps between transcripts of the same gene.
        Map<KeyGene, List<Region>> groups = new LinkedHashMap<>();

        for (Region r : regions) {
            if (r == null) {
                continue;
            }

            String gene = (r.gene == null) ? "" : r.gene.trim();
            if (gene.isEmpty()) {
                gene = "NA"; // чтобы не потерять записи без /gene
            }
            int a = Math.min(r.start, r.end);
            int b = Math.max(r.start, r.end);

            Region rn = new Region(r.transcriptId, gene, r.type, r.idx, a, b, r.strand);
            KeyGene k = new KeyGene(gene, rn.type, rn.strand);
            groups.computeIfAbsent(k, kk -> new ArrayList<>()).add(rn);
        }

        List<Region> out = new ArrayList<>();

        for (Map.Entry<KeyGene, List<Region>> e : groups.entrySet()) {
            List<Region> list = e.getValue();

            // Sort by coordinates: smallest -> largest
            list.sort(Comparator.<Region>comparingInt(x -> x.start).thenComparingInt(x -> x.end));

            int gap = mergeTouching ? 1 : 0;

            List<int[]> mergedSegs = new ArrayList<>();
            int curS = Integer.MIN_VALUE;
            int curE = Integer.MIN_VALUE;

            for (Region r : list) {
                if (mergedSegs.isEmpty()) {
                    curS = r.start;
                    curE = r.end;
                    mergedSegs.add(new int[]{curS, curE});
                    continue;
                }

                if (r.start <= curE + gap) {
                    curE = Math.max(curE, r.end);
                    mergedSegs.get(mergedSegs.size() - 1)[1] = curE;
                } else {
                    curS = r.start;
                    curE = r.end;
                    mergedSegs.add(new int[]{curS, curE});
                }
            }

            // Group metadata
            String gene = e.getKey().gene;
            String type = e.getKey().type;
            char strand = e.getKey().strand;

            // transcript_id loses its meaning after merging by gene
            String transcriptId = "MULTI";

            if ("EXON".equals(type)) {
                int idx = 1;
                for (int[] s : mergedSegs) {
                    out.add(new Region(transcriptId, gene, "EXON", String.valueOf(idx), s[0], s[1], strand));
                    idx++;
                }
            } else {
// PROMOTER: can become one or more after union (if different TSS are far apart)
                for (int[] s : mergedSegs) {
                    out.add(new Region(transcriptId, gene, "PROMOTER", "P", s[0], s[1], strand));
                }
            }
        }

// Global sorting of the final list by coordinates (as requested)
        out.sort(Comparator.<Region>comparingInt(r -> r.start).thenComparingInt(r -> r.end));

        return out;
    }

    static final class KeyGene {

        final String gene;
        final String type;
        final char strand;

        KeyGene(String gene, String type, char strand) {
            this.gene = (gene == null) ? "" : gene;
            this.type = (type == null) ? "" : type;
            this.strand = strand;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof KeyGene)) {
                return false;
            }
            KeyGene k = (KeyGene) o;
            return strand == k.strand && gene.equals(k.gene) && type.equals(k.type);
        }

        @Override
        public int hashCode() {
            int h = gene.hashCode();
            h = 31 * h + type.hashCode();
            h = 31 * h + (int) strand;
            return h;
        }
    }

    static List<Region> mergeAndDedup(List<Region> regions, boolean mergeTouching) {
        // Group by (transcript_id, type, strand).
        // If you want to group by gene instead of transcript_id, change the key.
        Map<Key, List<Region>> groups = new LinkedHashMap<>();
        for (Region r : regions) {
            if (r == null) {
                continue;
            }
            int a = Math.min(r.start, r.end);
            int b = Math.max(r.start, r.end);

            Region rn = new Region(r.transcriptId, r.gene, r.type, r.idx, a, b, r.strand);
            Key k = new Key(rn.transcriptId, rn.type, rn.strand);
            groups.computeIfAbsent(k, kk -> new ArrayList<>()).add(rn);
        }

        List<Region> out = new ArrayList<>();

        for (Map.Entry<Key, List<Region>> e : groups.entrySet()) {
            List<Region> list = e.getValue();

            // sort by start then end
            list.sort(Comparator.<Region>comparingInt(x -> x.start).thenComparingInt(x -> x.end));

            // merge + automatically delete duplicates
            List<Region> merged = new ArrayList<>();
            int curS = Integer.MIN_VALUE;
            int curE = Integer.MIN_VALUE;

            // To avoid losing metadata gene/transcript/type/strand — we take it from the first region in the group
            Region meta = list.get(0);

            for (Region r : list) {
                if (merged.isEmpty()) {
                    curS = r.start;
                    curE = r.end;
                    merged.add(r); // temporarily, we will replace it correctly below
                    continue;
                }

                int gap = mergeTouching ? 1 : 0;
                if (r.start <= curE + gap) {
                    // overlap (или touching)
                    curE = Math.max(curE, r.end);
                } else {
                    // record the previous interval
                    merged.set(merged.size() - 1, new Region(meta.transcriptId, meta.gene, meta.type, meta.idx, curS, curE, meta.strand));
                    // start new
                    curS = r.start;
                    curE = r.end;
                    merged.add(r);
                }
            }
            // finalise the last one
            merged.set(merged.size() - 1, new Region(meta.transcriptId, meta.gene, meta.type, meta.idx, curS, curE, meta.strand));

            // EXON reindexing after merge
            if ("EXON".equals(meta.type)) {
                int idx = 1;
                for (Region r : merged) {
                    out.add(new Region(r.transcriptId, r.gene, r.type, String.valueOf(idx), r.start, r.end, r.strand));
                    idx++;
                }
            } else {
                // PROMOTER: one (or several, if there were many ncRNAs) — idx can be left as "P"
                for (Region r : merged) {
                    String idxVal = (r.idx == null || r.idx.isBlank()) ? "P" : r.idx;
                    out.add(new Region(r.transcriptId, r.gene, r.type, idxVal, r.start, r.end, r.strand));
                }
            }
        }
        return out;
    }

    static final class Key {

        final String transcriptId;
        final String type;
        final char strand;

        Key(String transcriptId, String type, char strand) {
            this.transcriptId = (transcriptId == null) ? "" : transcriptId;
            this.type = (type == null) ? "" : type;
            this.strand = strand;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof Key)) {
                return false;
            }
            Key k = (Key) o;
            return strand == k.strand
                    && transcriptId.equals(k.transcriptId)
                    && type.equals(k.type);
        }

        @Override
        public int hashCode() {
            int h = transcriptId.hashCode();
            h = 31 * h + type.hashCode();
            h = 31 * h + (int) strand;
            return h;
        }
    }

    private static void processFeature(Feature f, int[] seqLenBox, int promoterLen, List<Region> out) throws IOException {
        // infer seqLen from source feature if LOCUS was not parsed
        if ("source".equals(f.key) && seqLenBox[0] < 0) {
            LocationParsed lp = parseLocation(f.location.toString());
            int maxEnd = lp.segments.stream().mapToInt(a -> a[1]).max().orElse(-1);
            if (maxEnd > 0) {
                seqLenBox[0] = maxEnd;
            }
            return;
        }

        if (!"ncRNA".equals(f.key)) {
            return;
        }

        String gene = f.getFirst("gene");
        String transcriptId = f.getFirst("transcript_id");

        LocationParsed lp = parseLocation(f.location.toString());
        List<int[]> segs = new ArrayList<>(lp.segments);
        if (segs.isEmpty()) {
            return;
        }

        // Sort in transcription order (useful for exon indexing)
        if (lp.strand == '+') {
            segs.sort(Comparator.comparingInt(a -> a[0]));
        } else {
            segs.sort((a, b) -> Integer.compare(b[1], a[1]));
        }

        // EXONS
        int exonIdx = 1;
        for (int[] s : segs) {
            out.add(new Region(transcriptId, gene, "EXON", String.valueOf(exonIdx),
                    s[0], s[1], lp.strand));
            exonIdx++;
        }

        // PROMOTER (computed)
        if (promoterLen > 0) {
            int seqLen = seqLenBox[0]; // may be -1 if unknown
            Region prom = computePromoter(transcriptId, gene, lp.strand, segs, promoterLen, seqLen);
            if (prom != null) {
                out.add(prom);
            }
        }
    }

    /**
     * Compute promoter coordinates upstream of TSS.
     */
    private static Region computePromoter(String transcriptId, String gene, char strand,
            List<int[]> segs, int promoterLen, int seqLen) {
        if (segs == null || segs.isEmpty()) {
            return null;
        }

        if (strand == '+') {
            int tss = segs.stream().mapToInt(a -> a[0]).min().orElse(-1);
            if (tss <= 1) {
                return null;
            }

            int start = Math.max(1, tss - promoterLen);
            int end = tss - 1;

            // clamp (seqLen not needed on + for upstream, but harmless)
            if (seqLen > 0) {
                start = Math.max(1, Math.min(start, seqLen));
                end = Math.max(1, Math.min(end, seqLen));
            }
            if (start <= end) {
                return new Region(transcriptId, gene, "PROMOTER", "P", start, end, strand);
            }
            return null;
        } else {
            int tss = segs.stream().mapToInt(a -> a[1]).max().orElse(-1);
            if (tss <= 0) {
                return null;
            }

            int start = tss + 1;
            int end = tss + promoterLen;

            // clamp to sequence end if known
            if (seqLen > 0) {
                start = Math.max(1, Math.min(start, seqLen));
                end = Math.max(1, Math.min(end, seqLen));
            }

            if (start <= end) {
                return new Region(transcriptId, gene, "PROMOTER", "P", start, end, strand);
            }
            return null;
        }
    }

    // ---------------------------- GenBank FEATURES parser ----------------------------
    private static class Feature {

        final String key;
        final StringBuilder location = new StringBuilder();
        final Map<String, List<String>> quals = new LinkedHashMap<>();

        boolean seenQualifier = false;
        String lastQualKey = null;
        StringBuilder lastQualValue = null;
        boolean lastQualOpenQuote = false;

        Feature(String key, String firstRest) {
            this.key = key;
            if (firstRest != null && !firstRest.isBlank()) {
                location.append(firstRest.trim());
            }
        }

        void consumeContinuation(String rest) {
            if (rest == null || rest.isBlank()) {
                return;
            }

            if (rest.startsWith("/")) {
                seenQualifier = true;
                parseQualifierLine(rest);
                return;
            }

            if (seenQualifier && lastQualKey != null && lastQualValue != null) {
                appendQualifierContinuation(rest);
                return;
            }

            if (!location.isEmpty()) {
                location.append(" ");
            }
            location.append(rest.trim());
        }

        String getFirst(String key) {
            List<String> v = quals.get(key);
            if (v == null || v.isEmpty()) {
                return "";
            }
            return v.get(0);
        }

        private void parseQualifierLine(String qline) {
            String q = qline.substring(1).trim();
            String k, v;

            int eq = q.indexOf('=');
            if (eq >= 0) {
                k = q.substring(0, eq).trim();
                v = q.substring(eq + 1).trim();
            } else {
                k = q.trim();
                v = "";
            }

            boolean openQuote = false;
            if (v.startsWith("\"")) {
                v = v.substring(1);
                if (v.endsWith("\"") && v.length() >= 1) {
                    v = v.substring(0, v.length() - 1);
                    openQuote = false;
                } else {
                    openQuote = true;
                }
            }

            quals.computeIfAbsent(k, kk -> new ArrayList<>()).add(v);

            lastQualKey = k;
            lastQualOpenQuote = openQuote;
            lastQualValue = new StringBuilder(v);

            // sync last stored value
            List<String> list = quals.get(k);
            list.set(list.size() - 1, lastQualValue.toString());
        }

        private void appendQualifierContinuation(String rest) {
            String s = rest.trim();

            if (lastQualOpenQuote) {
                if (s.endsWith("\"")) {
                    s = s.substring(0, s.length() - 1);
                    lastQualOpenQuote = false;
                }
            }

            if (!s.isBlank()) {
                if (lastQualValue.length() > 0) {
                    lastQualValue.append(" ");
                }
                lastQualValue.append(s);
            }

            List<String> list = quals.get(lastQualKey);
            if (list != null && !list.isEmpty()) {
                list.set(list.size() - 1, lastQualValue.toString());
            }
        }
    }

    // ---------------------------- Location parser ----------------------------
    private static class LocationParsed {

        final char strand;          // '+' or '-'
        final List<int[]> segments; // [start,end], start<=end

        LocationParsed(char strand, List<int[]> segments) {
            this.strand = strand;
            this.segments = segments;
        }
    }

    private static LocationParsed parseLocation(String loc) throws IOException {
        if (loc == null) {
            return new LocationParsed('+', List.of());
        }
        String s = loc.trim();

        char strand = '+';
        if (s.startsWith("complement(")) {
            strand = '-';
            s = insideParens(s);
        }
        s = s.trim();

        if (s.startsWith("join(") || s.startsWith("order(")) {
            String inner = insideParens(s);
            List<String> parts = splitTopLevelCommas(inner);
            List<int[]> segs = new ArrayList<>();
            for (String p : parts) {
                segs.addAll(parseSimpleOrNestedSegment(p.trim()));
            }
            return new LocationParsed(strand, segs);
        }

        return new LocationParsed(strand, parseSimpleOrNestedSegment(s));
    }

    private static List<int[]> parseSimpleOrNestedSegment(String token) throws IOException {
        String t = token.trim();

        if (t.startsWith("complement(")) {
            t = insideParens(t).trim();
        }

        t = t.replace("<", "").replace(">", "").replace(" ", "");

        if (t.startsWith("join(") || t.startsWith("order(")) {
            String inner = insideParens(t);
            List<String> parts = splitTopLevelCommas(inner);
            List<int[]> segs = new ArrayList<>();
            for (String p : parts) {
                segs.addAll(parseSimpleOrNestedSegment(p));
            }
            return segs;
        }

        int dd = t.indexOf("..");
        if (dd >= 0) {
            int a = Integer.parseInt(t.substring(0, dd));
            int b = Integer.parseInt(t.substring(dd + 2));
            return List.of(new int[]{Math.min(a, b), Math.max(a, b)});
        }

        int caret = t.indexOf('^');
        if (caret >= 0) {
            int a = Integer.parseInt(t.substring(0, caret));
            int b = Integer.parseInt(t.substring(caret + 1));
            return List.of(new int[]{Math.min(a, b), Math.max(a, b)});
        }

        int x = Integer.parseInt(t);
        return List.of(new int[]{x, x});
    }

    private static String insideParens(String s) throws IOException {
        int l = s.indexOf('(');
        int r = s.lastIndexOf(')');
        if (l < 0 || r < 0 || r <= l) {
            throw new IOException("Bad location: " + s);
        }
        return s.substring(l + 1, r);
    }

    private static List<String> splitTopLevelCommas(String s) {
        List<String> out = new ArrayList<>();
        int depth = 0;
        StringBuilder cur = new StringBuilder();
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '(') {
                depth++;
            }
            if (c == ')') {
                depth--;
            }
            if (c == ',' && depth == 0) {
                out.add(cur.toString());
                cur.setLength(0);
            } else {
                cur.append(c);
            }
        }
        if (cur.length() > 0) {
            out.add(cur.toString());
        }
        return out;
    }

    // ---------------------------- helpers ----------------------------
    private static String safeSubstring(String s, int from, int to) {
        if (s == null) {
            return "";
        }
        int a = Math.max(0, Math.min(from, s.length()));
        int b = Math.max(0, Math.min(to, s.length()));
        if (b < a) {
            return "";
        }
        return s.substring(a, b);
    }

    private static Integer parseLocusLength(String locusLine) {
        if (locusLine == null) {
            return null;
        }
        String[] tok = locusLine.trim().split("\\s+");
        for (int i = 0; i < tok.length - 1; i++) {
            if ("bp".equalsIgnoreCase(tok[i + 1])) {
                try {
                    return Integer.parseInt(tok[i]);
                } catch (NumberFormatException ignored) {
                }
            }
        }
        return null;
    }

    private static String nvl(String s) {
        return (s == null) ? "" : s;
    }

    // ---------------------------- output record ----------------------------
    public static class Region {

        public final String transcriptId;
        public final String gene;
        public final String type;   // EXON or PROMOTER
        public final String idx;    // 1..N for exons, "P" for promoter
        public final int start;
        public final int end;
        public final char strand;

        public Region(String transcriptId, String gene, String type, String idx, int start, int end, char strand) {
            this.transcriptId = transcriptId;
            this.gene = gene;
            this.type = type;
            this.idx = idx;
            this.start = start;
            this.end = end;
            this.strand = strand;
        }
    }
}
