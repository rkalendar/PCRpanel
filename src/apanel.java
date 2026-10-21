import java.nio.file.*;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class apanel {

    // ========================== CONFIG CLASS ==========================
    /**
     * Holds all primer-design parameters parsed from the configuration file.
     * Replaces 15+ loose local variables with a single, passable object.
     */
    static class PrimerConfig {

        int minPcr = 60;
        int maxPcr = 600;
        int minTm = 60;
        int maxTm = 62;
        int minLen = 18;
        int maxLen = 25;
        int prLap = 12;   // 0-18
        int minLc = 75;
        String e5 = "n";
        String e3 = "w";
        String fTail = "";
        String rTail = "";
        String outPath = null;
        public String[] tagMasks = {"*.*"};
        public String[] genomeMasks = {"*.*"};
        boolean homology = false;
        boolean multiplex = true;
    }

    // ========================== HELPERS ==========================
    /**
     * Clamp value into [min, max].
     */
    private static int clamp(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }

    /**
     * Extract the value portion after the FIRST '=' on the line. Returns "" if
     * no '=' is found.
     */
    private static String valueAfterEquals(String line) {
        int idx = line.indexOf('=');
        return (idx >= 0) ? line.substring(idx + 1).trim() : "";
    }

    /**
     * Parse an integer parameter from a config line that starts with the given
     * key. Returns null if the line doesn't contain the key.
     */
    private static Integer parseIntParam(String line, String key) {
        if (!line.contains(key)) {
            return null;
        }
        System.out.println(line);
        return StrToInt(valueAfterEquals(line));
    }

    /**
     * Build the concatenation of all sequences EXCEPT the one at index
     * {@code excludeIdx}.
     */
    private static String buildOtherSequences(List<String> sequences, int excludeIdx) {
        int totalLen = 0;
        for (int j = 0; j < sequences.size(); j++) {
            if (j == excludeIdx) {
                continue;
            }
            String s = sequences.get(j);
            if (s != null) {
                totalLen += s.length();
            }
        }
        StringBuilder sb = new StringBuilder(totalLen);
        for (int j = 0; j < sequences.size(); j++) {
            if (j == excludeIdx) {
                continue;
            }
            String s = sequences.get(j);
            if (s != null) {
                sb.append(s);
            }
        }
        return sb.toString();
    }

    /**
     * Read primers from a file, returning only those longer than 8 characters.
     */
    private static String[] readPrimers(String primerFile) {
        if (primerFile == null || primerFile.length() <= 1) {
            return new String[0];
        }
        List<String> result = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(primerFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim().toLowerCase();
                if (!line.isBlank()) {
                    String p = PrimerReturn(line);
                    if (p.length() > 8) {
                        result.add(p);
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Warning: could not read primer file '" + primerFile + "': " + e.getMessage());
        }
        return result.toArray(new String[0]);
    }

    /**
     * Read and concatenate reference FASTA file(s) into a single cleaned
     * sequence string.
     */
    private static String readReferenceSequences(List<String> refFiles) {
        if (refFiles.isEmpty()) {
            return "";
        }

        StringBuilder refBuilder = new StringBuilder(1 << 20);
        for (String refFile : refFiles) {
            try (BufferedReader br = new BufferedReader(
                    new InputStreamReader(
                            refFile.endsWith(".gz")
                            ? new java.util.zip.GZIPInputStream(new FileInputStream(refFile))
                            : new FileInputStream(refFile),
                            java.nio.charset.StandardCharsets.UTF_8), 1 << 16)) {

                String line;
                StringBuilder seqBuilder = new StringBuilder(1 << 20);
                while ((line = br.readLine()) != null) {
                    if (!line.isEmpty() && line.charAt(0) != '>') {
                        seqBuilder.append(line);
                    }
                }
                refBuilder.append(dna.DNA(seqBuilder.toString()));

            } catch (IOException e) {
                System.err.println("Warning: could not read reference file '" + refFile + "': " + e.getMessage());
            }
        }
        return refBuilder.toString();
    }

    /**
     * Load sequence + name + exons from a tagfile, with FASTA fallback. Returns
     * null if no valid sequence is found.
     */
    private static class SeqRecord {

        final String sequence;
        final String name;
        final int[] exons;
        final List<ExonPromoterExtractor.Region> regions; // annotation of each exon pair; null for FASTA

        SeqRecord(String sequence, String name, int[] exons, List<ExonPromoterExtractor.Region> regions) {
            this.sequence = sequence;
            this.name = name;
            this.exons = exons;
            this.regions = regions;
        }
    }

    private static SeqRecord loadSequenceRecord(String tagFile) throws Exception {
        ExonPromoterExtractor ex = new ExonPromoterExtractor(tagFile);
        String sequence = ex.getSequence();
        String name = ex.getRecordId();

        if (sequence != null && !sequence.isBlank()) {
            List<ExonPromoterExtractor.Region> regs = ex.getRegions();
            int[] exons = regs.stream()
                    .flatMapToInt(r -> java.util.stream.IntStream.of(r.start, r.end))
                    .toArray();
            return new SeqRecord(sequence, name, exons, regs);
        }

        // FASTA fallback — use last valid record
        List<FastaIO.FastaRecord> records
                = FastaIO.readAllToList(java.util.Collections.singletonList(tagFile));
        for (FastaIO.FastaRecord rec : records) {
            if (rec != null && rec.sequence != null && !rec.sequence.isBlank()) {
                String recName = (rec.id != null && !rec.id.isBlank()) ? rec.id : name;
                int[] exons = new int[]{0, rec.sequence.length()};
                // Note: returns LAST valid record (preserving original behavior)
                sequence = rec.sequence;
                name = recName;
                return new SeqRecord(sequence, name, exons, null);
            }
        }
        return null; // nothing found
    }

    private static String[] splitMasks(String s) {
        if (s == null || s.isBlank()) {
            return new String[]{"*.*"};
        }
        return java.util.Arrays.stream(s.split("[,;|]"))
                .map(String::trim)
                .filter(p -> !p.isEmpty())
                // "fna" -> "*.fna"; "*.fna" / "name_*.fa" оставляем как есть
                .map(p -> p.contains("*") || p.contains("?") ? p : "*." + p)
                .toArray(String[]::new);
    }

    private static void dedupInPlace(List<String> list) {
        java.util.LinkedHashSet<String> set = new java.util.LinkedHashSet<>(list);
        list.clear();
        list.addAll(set);
    }

    // ========================== CONFIG PARSING ==========================
    /**
     * Parse the configuration file and populate all settings.
     */
    private static PrimerConfig parseConfig(String infile,
            List<String> tagfiles,
            List<String> reffiles,
            List<String> genomefiles,
            /* out */ String[] primerFileHolder) throws Exception {
        PrimerConfig cfg = new PrimerConfig();

        // NEW: собираем папки во временные списки — раскроем после того, как узнаем маски
        List<String> tagFolders = new ArrayList<>();
        List<String> genomeFolders = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String cline = line;
                line = line.toLowerCase();

                // --- Path-based parameters ---
                Path folderPath = validatePath(cline, "folder_path=", true, false);
                if (folderPath != null) {
                    String s = cline.substring(cline.indexOf("folder_path=") + "folder_path=".length()).trim();
                    System.out.println("Target_Folder_path=" + s);
                    tagFolders.add(s);                                  // CHANGED: было BatchReadFiles здесь
                }

                Path genomePath = validatePath(cline, "genome_path=", true, false);
                if (genomePath != null) {
                    String s = cline.substring(cline.indexOf("genome_path=") + "genome_path=".length()).trim();
                    System.out.println("Genome_path=" + genomePath);
                    genomeFolders.add(s);                               // CHANGED
                }

                Path outDir = validatePath(cline, "folder_out=", true, true);
                if (outDir != null) {
                    cfg.outPath = cline.substring(cline.indexOf("folder_out=") + "folder_out=".length()).trim();
                    cfg.outPath = stripQuotes(cfg.outPath);
                    outDir = OutDirUtil.prepareOutputDir(cfg.outPath);
                    System.out.println("Output_Folder_path=" + outDir);
                    cfg.outPath = outDir.toString();
                }

                Path inputFile = validatePath(cline, "target_path=", false, false);
                if (inputFile != null) {
                    System.out.println("Input_file= " + inputFile);
                    tagfiles.add(inputFile.toString());
                }

                Path refFile = validatePath(cline, "reference_path=", false, false);
                if (refFile != null) {
                    System.out.println("Reference_file_path= " + refFile);
                    reffiles.add(refFile.toString());
                }

                // NEW: маски файлов
                if (line.startsWith("target_mask=") || line.startsWith("folder_mask=")) {
                    cfg.tagMasks = splitMasks(valueAfterEquals(line));
                    System.out.println("Target_mask=" + String.join(",", cfg.tagMasks));
                }
                if (line.startsWith("genome_mask=")) {
                    cfg.genomeMasks = splitMasks(valueAfterEquals(line));
                    System.out.println("Genome_mask=" + String.join(",", cfg.genomeMasks));
                }
                // Опционально: общая маска, если пользователь хочет одну на всё
                if (line.startsWith("file_mask=")) {
                    String[] m = splitMasks(valueAfterEquals(line));
                    cfg.tagMasks = m;
                    cfg.genomeMasks = m;
                    System.out.println("File_mask=" + String.join(",", m));
                }

                // --- Boolean flags ---
                if (line.contains("homology=true")) {
                    cfg.homology = true;
                    System.out.println("Designing common primers only based on shared sequences between different files.");
                }
                if (line.contains("multiplex=")) {
                    cfg.multiplex = !line.contains("multiplex=false");
                    System.out.println(cfg.multiplex ? "Multiplex two panel design" : "Single-plex panel design");
                }

                // --- String parameters ---
                if (line.startsWith("3end=")) {
                    cfg.e3 = valueAfterEquals(line);
                    System.out.println("3End=" + cfg.e3);
                }
                if (line.startsWith("5end=")) {
                    cfg.e5 = valueAfterEquals(line);
                    System.out.println("5End=" + cfg.e5);
                }
                if (line.startsWith("target_primers=")) {
                    primerFileHolder[0] = valueAfterEquals(line);
                    System.out.println("Target_primers=" + primerFileHolder[0]);
                }
                if (line.startsWith("forwardtail=")) {
                    cfg.fTail = dna.DNA(valueAfterEquals(line));
                    System.out.println("ForwardTail=" + cfg.fTail);
                }
                if (line.startsWith("reversetail=")) {
                    cfg.rTail = dna.DNA(valueAfterEquals(line));
                    System.out.println("ReverseTail=" + cfg.rTail);
                }

                // --- Numeric parameters with clamping ---
                Integer v;
                if ((v = parseIntParam(line, "minpcr=")) != null) {
                    cfg.minPcr = clamp(v, 30, 5000);
                }
                if ((v = parseIntParam(line, "maxpcr=")) != null) {
                    cfg.maxPcr = clamp(v, 30, 50000);
                }
                if ((v = parseIntParam(line, "minlen=")) != null) {
                    cfg.minLen = clamp(v, 12, 80);
                }
                if ((v = parseIntParam(line, "maxlen=")) != null) {
                    cfg.maxLen = clamp(v, 12, 100);
                }
                if ((v = parseIntParam(line, "mintm=")) != null) {
                    cfg.minTm = clamp(v, 40, 75);
                }
                if ((v = parseIntParam(line, "maxtm=")) != null) {
                    cfg.maxTm = clamp(v, 40, 80);
                }
                if ((v = parseIntParam(line, "minlc=")) != null) {
                    cfg.minLc = clamp(v, 0, 100);
                    System.out.println("minLC=" + cfg.minLc);
                }
                if ((v = parseIntParam(line, "prlap=")) != null) {
                    cfg.prLap = clamp(v, 0, 18);
                    System.out.println("prLap=" + cfg.prLap);
                }
                if ((v = parseIntParam(line, "maskedbases=")) != null) {
                    PrimerDesign.MASK_TOLERANCE = clamp(v, 0, 40);
                    System.out.println("maskedBases=" + PrimerDesign.MASK_TOLERANCE);
                }
                if ((v = parseIntParam(line, "repeatkmer=")) != null) {
                    PrimerDesign.GENOME_KMER_SIZE = clamp(v, 8, 64);
                    System.out.println("repeatKmer=" + PrimerDesign.GENOME_KMER_SIZE);
                }
                if ((v = parseIntParam(line, "repeatmincount=")) != null) {
                    MaskingSequence.MIN_COUNT = clamp(v, 2, 1000);
                    System.out.println("repeatMinCount=" + MaskingSequence.MIN_COUNT);
                }
            }
        } catch (IOException e) {
            System.err.println("Error reading config file '" + infile + "': " + e.getMessage());
        }

        // NEW: раскрываем папки уже с учётом масок
        BatchReadFiles brf = new BatchReadFiles();
        for (String folder : tagFolders) {
            for (String mask : cfg.tagMasks) {
                String[] files = brf.BatchReadFiles(new String[]{folder, mask});
                if (files != null && files.length > 0) {
                    tagfiles.addAll(List.of(files));
                    System.out.println("Matched " + files.length + " file(s) for " + folder + " / " + mask);
                }
            }
        }
        for (String folder : genomeFolders) {
            for (String mask : cfg.genomeMasks) {
                String[] files = brf.BatchReadFiles(new String[]{folder, mask});
                if (files != null && files.length > 0) {
                    genomefiles.addAll(List.of(files));
                    System.out.println("Matched " + files.length + " file(s) for " + folder + " / " + mask);
                }
            }
        }
        dedupInPlace(tagfiles);
        dedupInPlace(genomefiles);

        return cfg;
    }

    // ========================== SINGLE-PLEX ==========================
    private static void runSingleplexAll(LongIntHashMap map, List<String> tagfiles,
            PrimerConfig cfg, String refsequence,
            String[] primers) throws Exception {
        final String ref = "";
        for (String tagfile : tagfiles) {
            ExonPromoterExtractor ex = new ExonPromoterExtractor(tagfile);
            String sequence = ex.getSequence();
            String name = ex.getRecordId();

            // FASTA fallback
            if (sequence == null || sequence.isBlank()) {
                List<FastaIO.FastaRecord> records
                        = FastaIO.readAllToList(java.util.Collections.singletonList(tagfile));
                if (records.isEmpty()) {
                    continue;
                }

                for (FastaIO.FastaRecord rec : records) {
                    if (rec == null || rec.sequence == null || rec.sequence.isBlank()) {
                        continue;
                    }
                    String recName = (rec.id != null && !rec.id.isBlank()) ? rec.id : name;
                    int[] exons = new int[]{0, rec.sequence.length()};
                    RunSingleplexPanelDesign(map, tagfile, cfg.outPath, rec.sequence, recName,
                            refsequence, ref, cfg.homology, primers, exons,
                            cfg.minPcr, cfg.maxPcr, cfg.minLen, cfg.maxLen,
                            cfg.minTm, cfg.maxTm, cfg.minLc, cfg.prLap,
                            cfg.fTail, cfg.rTail, cfg.e5, cfg.e3, null);
                }
                continue;
            }

            // GenBank / extracted regions path
            List<ExonPromoterExtractor.Region> regs = ex.getRegions();
            int[] exons = toExonArray(regs, sequence.length());
            RunSingleplexPanelDesign(map, tagfile, cfg.outPath, sequence, name,
                    refsequence, ref, cfg.homology, primers, exons,
                    cfg.minPcr, cfg.maxPcr, cfg.minLen, cfg.maxLen,
                    cfg.minTm, cfg.maxTm, cfg.minLc, cfg.prLap,
                    cfg.fTail, cfg.rTail, cfg.e5, cfg.e3, regs);
        }
    }

    // ========================== MULTIPLEX ==========================
    private static void runMultiplexAll(LongIntHashMap map, List<String> tagfiles,
            PrimerConfig cfg, String refsequence,
            String[] primers) throws Exception {
        List<String> sequences = new ArrayList<>();
        List<String> nameSequences = new ArrayList<>();
        List<int[]> exonsSequences = new ArrayList<>();
        List<String> fileSequences = new ArrayList<>();
        List<List<ExonPromoterExtractor.Region>> regionsSequences = new ArrayList<>();

        // Phase 1: load all sequences
        for (String tagfile : tagfiles) {
            SeqRecord rec = loadSequenceRecord(tagfile);
            if (rec == null) {
                continue;
            }
            sequences.add(rec.sequence);
            nameSequences.add(rec.name);
            exonsSequences.add(rec.exons);
            fileSequences.add(tagfile);
            regionsSequences.add(rec.regions);
        }

        // Phase 2: design primers file-by-file
        List<MultiplexTarget> designs = new ArrayList<>();
        int nseq = sequences.size();
        for (int i = 0; i < nseq; i++) {
            String sequence = sequences.get(i);
            if (sequence == null || sequence.isEmpty()) {
                continue;
            }

            String name = nameSequences.get(i);
            int[] exons = exonsSequences.get(i);
            String tagfile = fileSequences.get(i);
            String otherSeqs = buildOtherSequences(sequences, i);

            MultiplexTarget t = designMultiplexTarget(map, tagfile,
                    sequence, name, refsequence, otherSeqs,
                    cfg.homology, primers, exons,
                    cfg.minPcr, cfg.maxPcr, cfg.minLen, cfg.maxLen,
                    cfg.minTm, cfg.maxTm, cfg.minLc, cfg.prLap,
                    cfg.fTail, cfg.rTail, cfg.e5, cfg.e3);
            t.regions = regionsSequences.get(i);
            designs.add(t);
        }

        // Phase 3: one panel for all files, which share the two tubes. Assembling
        // them together lets the regions with the fewest primer pairs, in any file,
        // be placed first; designing each file against the primers of the files
        // before it left the last files without primers for some of their regions.
        List<Pair> pcrcol = new ArrayList<>();
        List<Pair> alternatives = new ArrayList<>();
        List<int[]> targets = new ArrayList<>();
        List<RegionNames> names = new ArrayList<>();
        for (int s = 0; s < designs.size(); s++) {
            MultiplexTarget t = designs.get(s);
            for (Pair p : t.pcrcol) {
                p.seq = s;
            }
            for (Pair p : t.allPairs) {
                p.seq = s;
            }
            pcrcol.addAll(t.pcrcol);
            alternatives.addAll(t.allPairs);
            targets.add(t.cexons);
            names.add(new RegionNames(t.exons, t.regions));
        }
        PanelsCollector panels = new PanelsCollector(pcrcol, alternatives);
        panels.setRegionNames(names);
        panels.CombinePairsPanel(targets, 2);

        // Phase 4: reports, file by file
        for (int s = 0; s < designs.size(); s++) {
            reportMultiplexTarget(designs.get(s), s, panels, cfg.outPath, cfg.maxPcr);
        }
    }

    // ========================== HELP / USAGE ==========================

    /** Application name used in help/usage output. */
    private static final String APP_NAME = "PCRpanel";

    /**
     * Returns {@code true} if the supplied command-line argument is a request
     * for help, e.g. {@code "/?"}, {@code "-help"}, {@code "--help"},
     * {@code "-h"}, {@code "?"} or {@code "help"} (case-insensitive).
     */
    private static boolean isHelpRequested(String arg) {
        if (arg == null) {
            return false;
        }
        switch (arg.trim().toLowerCase()) {
            case "/?":
            case "?":
            case "-?":
            case "-h":
            case "/h":
            case "-help":
            case "--help":
            case "/help":
            case "help":
                return true;
            default:
                return false;
        }
    }

    /** Format a single "key — description" line for the help screen. */
    private static String helpRow(String key, String desc) {
        return String.format("  %-24s %s%n", key, desc);
    }

    /**
     * Print program usage together with the full list of control commands that
     * can be placed in the configuration file. Triggered by {@code "/?"},
     * {@code "-help"} (and common variants), or when no configuration file is
     * supplied on the command line.
     */
    private static void printHelp() {
        String nl = System.lineSeparator();
        StringBuilder h = new StringBuilder(4096);

        h.append(APP_NAME).append(" - Custom Amplicon Panel Designer").append(nl);
        h.append("Designs custom multiplex / singleplex PCR amplicon panels for NGS and ONT.").append(nl);
        h.append("Documentation: https://primerdigital.com/tools/panel.html").append(nl);
        h.append(nl);

        h.append("USAGE").append(nl);
        h.append("  java -jar PCRpanel.jar <config-file>").append(nl);
        h.append("  java -jar PCRpanel.jar (-help | /?)").append(nl);
        h.append(nl);
        h.append("  <config-file>   Plain-text configuration file holding all control").append(nl);
        h.append("                  commands. Its name and extension do not matter").append(nl);
        h.append("                  (e.g. config.file, my_panel.conf, run.txt).").append(nl);
        h.append(nl);
        h.append("  For large reference genomes raise the JVM heap, e.g.:").append(nl);
        h.append("    java -Xms8g -Xmx32g -jar PCRpanel.jar config.file").append(nl);
        h.append(nl);

        h.append("CONFIGURATION FILE FORMAT").append(nl);
        h.append("  Write one command per line as:  key=value").append(nl);
        h.append("  Commands may appear in any order. Blank lines, lines starting with '#',").append(nl);
        h.append("  and any unrecognised text are ignored, so notes can be kept in the file.").append(nl);
        h.append("  Path commands (target_path, reference_path) may be repeated.").append(nl);
        h.append("  Keys are case-insensitive; out-of-range numbers are clamped to the limit.").append(nl);
        h.append(nl);

        h.append("INPUT / OUTPUT").append(nl);
        h.append(helpRow("target_path=FILE",    "Target GenBank/FASTA file. Repeatable."));
        h.append(helpRow("target_primers=FILE", "Existing primer/probe list to incorporate."));
        h.append(helpRow("reference_path=FILE", "Reference sequence file. Repeatable."));
        h.append(helpRow("folder_path=DIR",     "Folder of target files (subfolders included)."));
        h.append(helpRow("genome_path=DIR",     "Folder of reference-genome FASTA files (subfolders)."));
        h.append(helpRow("folder_out=DIR",      "Output folder. Created if missing; if it already"));
        h.append(helpRow("",                    "exists its CONTENTS ARE DELETED and replaced."));
        h.append(nl);

        h.append("FILE MASKS  (restrict which files are read from a folder)").append(nl);
        h.append(helpRow("target_mask=MASK",    "Mask for folder_path. Alias: folder_mask. Default *.*"));
        h.append(helpRow("genome_mask=MASK",    "Mask for genome_path. Default *.*"));
        h.append(helpRow("file_mask=MASK",      "Apply one mask to both folders."));
        h.append("    Forms: *.fna  |  fna (=> *.fna)  |  fna,fasta,gb  |  *.fna;*.gb").append(nl);
        h.append(nl);

        h.append("PANEL MODE").append(nl);
        h.append(helpRow("multiplex=true|false", "Two overlapping multiplex pools. Default true."));
        h.append(helpRow("",                     "Amplicons that conflict with both pools go to extra"));
        h.append(helpRow("",                     "tubes (Panel3, Panel4, ...), each free of conflicts,"));
        h.append(helpRow("",                     "instead of being dropped."));
        h.append(helpRow("homology=true|false",  "Consensus primers from regions shared across inputs."));
        h.append(helpRow("",                     "Default false."));
        h.append(nl);

        h.append("AMPLICON SIZE (bp)").append(nl);
        h.append(helpRow("minPCR=INT",          "Min amplicon size. Default 60.  Range 30-5000."));
        h.append(helpRow("maxPCR=INT",          "Max amplicon size. Default 600. Range 30-50000."));
        h.append(nl);

        h.append("PRIMER CONSTRAINTS").append(nl);
        h.append(helpRow("minLen=INT",          "Min primer length (nt). Default 18. Range 12-80."));
        h.append(helpRow("maxLen=INT",          "Max primer length (nt). Default 25. Range 12-100."));
        h.append(helpRow("minTm=INT",           "Min melting temperature (\u00B0C). Default 60. Range 40-75."));
        h.append(helpRow("maxTm=INT",           "Max melting temperature (\u00B0C). Default 62. Range 40-80."));
        h.append(helpRow("3end=CODE",           "3' terminal base (IUPAC). Default w."));
        h.append(helpRow("5end=CODE",           "5' terminal base (IUPAC). Default n (any base)."));
        h.append("    IUPAC: w=A/T  s=G/C  n=any;  r y m k b d h v also accepted.").append(nl);
        h.append(nl);

        h.append("ADAPTER TAILS").append(nl);
        h.append(helpRow("forwardtail=SEQ",     "5' adapter/tail prepended to forward primers."));
        h.append(helpRow("reversetail=SEQ",     "5' adapter/tail prepended to reverse primers."));
        h.append(nl);

        h.append("EXAMPLE (minimal configuration file)").append(nl);
        h.append("  target_path=/data/genes/NG_013019.gb").append(nl);
        h.append("  target_path=/data/genes/NG_011731.gb").append(nl);
        h.append("  multiplex=true").append(nl);
        h.append("  minPCR=250").append(nl);
        h.append("  maxPCR=500").append(nl);
        h.append("  minLen=18").append(nl);
        h.append("  maxLen=24").append(nl);
        h.append("  minTm=60").append(nl);
        h.append("  maxTm=62").append(nl);
        h.append("  3end=w").append(nl);
        h.append("  forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT").append(nl);
        h.append("  reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT").append(nl);
        h.append(nl);
        h.append("Then run:  java -jar PCRpanel.jar <your-config-file>").append(nl);
        h.append("See the README for full details, more examples, and best practices.").append(nl);

        System.out.println(h);
    }

    // ========================== MAIN ==========================
    public static void main(String[] args) throws IOException, Exception {
        // Show help for "/?", "-help", "--help", "-h" (and variants), or when
        // no configuration file is supplied on the command line.
        if (args.length == 0 || isHelpRequested(args[0])) {
            printHelp();
            return;
        }

        String infile = args[0];
        System.out.println("Current Directory: " + System.getProperty("user.dir"));
        System.out.println("Command-line arguments: " + infile);

        // Mutable lists populated by parseConfig
        List<String> tagfiles = new ArrayList<>();
        List<String> reffiles = new ArrayList<>();
        List<String> genomefiles = new ArrayList<>();
        String[] primerFileHolder = {""};  // single-element array to allow mutation

        // 1. Parse configuration
        PrimerConfig cfg = parseConfig(infile, tagfiles, reffiles, genomefiles, primerFileHolder);

        // 2. Read primers
        String[] primers = readPrimers(primerFileHolder[0]);

        // 3. Read reference sequences
        String refsequence = readReferenceSequences(reffiles);

        // 4. Build genome hash-map
        long startTime = System.nanoTime();
        LongIntHashMap map = GenomeReference(genomefiles);
        long duration = (System.nanoTime() - startTime) / 1_000_000_000;
        if (!map.isEmpty()) {
            System.out.println("Genome repeat searching, time taken: " + duration + " seconds\n");
            startTime = System.nanoTime();
        }

        // 5. Run design
        if (!cfg.multiplex) {
            runSingleplexAll(map, tagfiles, cfg, refsequence, primers);
        } else {
            runMultiplexAll(map, tagfiles, cfg, refsequence, primers);
        }

        duration = (System.nanoTime() - startTime) / 1_000_000_000;
        System.out.println("Time taken: " + duration + " seconds\n");
    }

    /**
     * Converts Region list to [start1,end1,start2,end2,...] with basic
     * sanitation.
     */
    private static int[] toExonArray(List<ExonPromoterExtractor.Region> regs, int seqLen) {
        if (regs == null || regs.isEmpty()) {
            return new int[]{0, Math.max(0, seqLen)};
        }

        int[] exons = new int[regs.size() * 2];
        int i = 0;

        for (ExonPromoterExtractor.Region r : regs) {
            int s = r.start;
            int e = r.end;

            // normalize order
            if (e < s) {
                int tmp = s;
                s = e;
                e = tmp;
            }

            // clamp to sequence bounds
            s = Math.max(0, Math.min(s, seqLen));
            e = Math.max(0, Math.min(e, seqLen));

            exons[i++] = s;
            exons[i++] = e;
        }

        return exons;
    }

    /**
     * One target file of a multiplex run, between its primer design and its
     * report. It keeps the design results rather than the PrimerDesign, which
     * holds the other files' sequences.
     */
    private static class MultiplexTarget {

        String tagfile;
        String seq;
        String name;      // record ID, as in the primer names
        int[] exons;      // regions as merged by the design
        int[] cexons;     // regions before merging
        List<ExonPromoterExtractor.Region> regions; // annotation of the regions; null for FASTA
        PrimersCollector[] fPrimersList;
        PrimersCollector[] rPrimersList;
        int[] widened;
        List<Pair> pcrcol;
        List<Pair> allPairs;
        StringBuilder sr; // primer list report, its header written at design time
    }

    private static MultiplexTarget designMultiplexTarget(LongIntHashMap map, String tagfile, String seq, String name, String refsequence, String fastaseq, Boolean homology, String[] listprimers, int[] exons, int minpcr, int maxpcr, int minlen, int maxlen, int mintm, int maxtm, int minlc, int prlap, String ftail, String rtail, String e5, String e3) {
        MultiplexTarget t = new MultiplexTarget();
        t.tagfile = tagfile;
        t.seq = seq;
        t.name = name;
        t.exons = exons;
        t.cexons = exons.clone();

        /*   ExonsDesigner pex = new ExonsDesigner(minpcr, maxpcr);
        List<ExonsDesigner.PcrFragment> fr = pex.buildPcrFragments(cexons);
        fr.forEach(System.out::println);
         */
        StringBuilder sr = new StringBuilder(100000);
        t.sr = sr;

        long startTime = System.nanoTime();
        System.out.println("Running...");
        System.out.println("\nTarget file name: " + tagfile);

        sr.append("Target file name: ").append(tagfile).append("\n");
        if (homology) {
            sr.append("Designing common primers only based on shared sequences between different files.\n");
        }
        sr.append("min PCR size=").append(minpcr).append("\n");
        sr.append("max PCR size=").append(maxpcr).append("\n");
        sr.append("min Tm=").append(mintm).append("\n");
        sr.append("max Tm=").append(maxtm).append("\n");
        sr.append("min length=").append(minlen).append("\n");
        sr.append("max length=").append(maxlen).append("\n");

        if (e3.length() > 0) {
            sr.append("3-end=").append(e3).append("\n");
        }
        if (e5.length() > 0) {
            sr.append("5-end=").append(e5).append("\n");
        }
        if (ftail.length() > 0) {
            sr.append("Forward tail=").append(ftail).append("\n");
        }
        if (rtail.length() > 0) {
            sr.append("Reverse tail=").append(rtail).append("\n");
        }

        PrimerDesign pd = new PrimerDesign(map, seq, name, minlen, maxlen, mintm, maxtm, minlc, exons, listprimers, refsequence, fastaseq, homology);
        t.pcrcol = pd.runDesign(ftail, rtail, e5, e3, prlap, minpcr, maxpcr);
        t.allPairs = pd.getAllPairs();
        t.fPrimersList = pd.getForwardPrimers();
        t.rPrimersList = pd.getReversePrimers();
        t.widened = pd.getWidenedFlanks();

        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Time taken: " + duration + " seconds\n\n");
        return t;
    }

    /**
     * Writes the reports of target {@code s}, whose amplicons the shared
     * {@code panels} carry.
     *
     * @param maxpcr maximum amplicon size, also the flank written around the
     *               regions left uncovered
     */
    private static void reportMultiplexTarget(MultiplexTarget t, int s, PanelsCollector panels, String outpath, int maxpcr) {
        String tagfile = t.tagfile;
        String seq = t.seq;
        int[] exons = t.exons;
        int[] cexons = t.cexons;
        List<Pair> pcrcol = t.pcrcol;

        Path in = Paths.get(tagfile);
        // String stem = fileStem(in).replaceFirst("\\.[0-9]+$", "");  // NG_029916.1
        String stem = fileStem(in);
        Path outDir = (outpath == null || outpath.isBlank()) ? in.getParent() : Paths.get(outpath);
        if (outDir == null) {
            outDir = Paths.get(".");
        }
        // Path outDir = base.resolve(stem);

        Path primerlistfile = outDir.resolve(stem + "_primers.txt");
        Path pcrcolfile = outDir.resolve(stem + "_panels.txt");
        Path panelprimerlistfile = outDir.resolve(stem + "_panelprimers.txt");

        try {
            StringBuilder sr = t.sr;
            StringBuilder sr1 = new StringBuilder(100000);
            StringBuilder srpanel = new StringBuilder(100000);

            System.out.println("\nTarget file name: " + tagfile);
            RegionNames names = new RegionNames(exons, t.regions);

            // Panel1 and Panel2 are the two multiplex tubes; Panel3 onwards, opened only
            // when needed, carry the amplicons that conflict with both
            List<List<Pair>> pcrpanels = new ArrayList<>();
            for (int k = 1; k <= panels.getPanelCount(); k++) {
                pcrpanels.add(poolOf(panels.getPanel(k), s));
            }
            for (int k = 0; k < pcrpanels.size(); k++) {
                appendPanel(sr1, srpanel, "Panel" + (k + 1), pcrpanels.get(k), seq, names);
            }
            String notes = panels.renderNotes(s);
            if (!notes.isEmpty()) {
                sr1.append(notes);
                System.out.println(notes);
            }
            CoverageReport cov = new CoverageReport(cexons, t.regions);
            for (int k = 0; k < pcrpanels.size(); k++) {
                if (!pcrpanels.get(k).isEmpty()) {
                    List<int[]> iv = new ArrayList<>();
                    for (Pair pcrx : pcrpanels.get(k)) {
                        iv.add(CoverageReport.interval(pcrx.fx5, 1 + pcrx.rx5));
                    }
                    cov.addPool("Panel" + (k + 1), iv);
                }
            }
            String covtxt = cov.render();
            if (!covtxt.isEmpty()) {
                sr1.append(covtxt);
                System.out.println(covtxt);
            }
            /*
            // pcrcol - here all pairs list with exons ID
            if (!pcrcol.isEmpty()) {
                sr1.append("List all amplicons:\n");
                for (Pair pcrx : pcrcol) {
                    sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
                    sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.rlc).append("\n");
                    sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp\n\n");
                }
            }
             */
            PrimersCollector[] fPrimersList = t.fPrimersList;
            PrimersCollector[] rPrimersList = t.rPrimersList;
            int[] widened = t.widened;
            int h = -1;
            sr.append("\nName\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)");
            for (int j = 0; j < exons.length - 1; j += 2) {
                h++;
                String region = regionLine(names, exons, cexons, widened, h);
                sr.append("\n").append(region).append("\n");
                System.out.println("\n" + region);

                PrimersCollector fPrimersList1 = fPrimersList[h];
                if (fPrimersList1.Amount() > 0) {
                    StringBuilder sr2 = new StringBuilder(100000);
                    double[] Tm = fPrimersList1.getTms();
                    double[] CG = fPrimersList1.getCGs();
                    int[] ln = fPrimersList1.getPrimerLengths();          // length
                    int[] lc = fPrimersList1.getPrimerLC();               // Linguistic_Complexity
                    // int[] x1 = fPrimersList1.getPrimerLocations();        // location x1
                    String[] primer = fPrimersList1.getPrimer();          // primer sequence
                    String[] primername = fPrimersList1.getpPrimerName(); // name

                    for (int i = 0; i < fPrimersList1.Amount(); i++) {
                        sr2.append(primername[i]).append("\t").append(primer[i]).append("\t").append(ln[i]).append("\t").append(String.format("%.1f", Tm[i])).append("\t").append(String.format("%.1f", CG[i])).append("\t").append(lc[i]).append("\n");
                    }
                    System.out.println(sr2);
                    sr.append(sr2);
                }

                PrimersCollector rPrimersList1 = rPrimersList[h];
                if (rPrimersList1.Amount() > 0) {
                    StringBuilder sr2 = new StringBuilder(100000);
                    double[] Tm = rPrimersList1.getTms();
                    double[] CG = rPrimersList1.getCGs();
                    int[] ln = rPrimersList1.getPrimerLengths();          // length
                    int[] lc = rPrimersList1.getPrimerLC();               // Linguistic_Complexity
                    //  int[] x1 = rPrimersList1.getPrimerLocations();        // location x1
                    String[] primer = rPrimersList1.getPrimer();          // primer sequence
                    String[] primername = rPrimersList1.getpPrimerName(); //name
                    for (int i = 0; i < rPrimersList1.Amount(); i++) {
                        sr2.append(primername[i]).append("\t").append(primer[i]).append("\t").append(ln[i]).append("\t").append(String.format("%.1f", Tm[i])).append("\t").append(String.format("%.1f", CG[i])).append("\t").append(lc[i]).append("\n");
                    }
                    System.out.println(sr2);
                    sr.append(sr2);
                }
            }

            try (FileWriter fileWriter = new FileWriter(primerlistfile.toFile())) {
                System.out.println("Saving the primer list report to a file: " + primerlistfile);
                fileWriter.write(sr.toString());
            }

            if (srpanel.length() > 0) {
                try (FileWriter fileWriter = new FileWriter(panelprimerlistfile.toFile())) {
                    System.out.println("Saving panels PCR primers combinations report to a file: " + panelprimerlistfile);
                    fileWriter.write(srpanel.toString());
                }
            }

            PanelPlot plot = new PanelPlot(t.name, seq.length(), cexons, exons, t.regions, t.widened, pcrpanels, true);
            if (!pcrcol.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(pcrcolfile.toFile())) {
                    System.out.println("Saving PCR primers combinations report to a file: " + pcrcolfile);
                    fileWriter.write(sr1.toString());
                }
                writePlot(outDir, stem, plot);
            }
            writeUncovered(outDir, stem, plot.uncoveredFasta(seq, maxpcr));

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    /**
     * Writes the parts of the targets left uncovered, with their flanks, as
     * FASTA; nothing when every target is covered in full.
     */
    private static void writeUncovered(Path outDir, String stem, String fasta) throws IOException {
        if (fasta.isEmpty()) {
            return;
        }
        Path file = outDir.resolve(stem + "_uncovered.fasta");
        System.out.println("Saving the regions not covered in full (FASTA) to a file: " + file);
        Files.writeString(file, fasta);
    }

    /** Writes the panel map (SVG) and the BED tracks of one target. */
    private static void writePlot(Path outDir, String stem, PanelPlot plot) throws IOException {
        Path svg = outDir.resolve(stem + "_panel.svg");
        System.out.println("Saving the panel map to a file: " + svg);
        Files.writeString(svg, plot.svg());
        Path bed = outDir.resolve(stem + "_panelprimers.bed");
        System.out.println("Saving the panel primers (BED) to a file: " + bed);
        Files.writeString(bed, plot.primersBed());
        bed = outDir.resolve(stem + "_amplicons.bed");
        System.out.println("Saving the amplicons (BED) to a file: " + bed);
        Files.writeString(bed, plot.ampliconsBed());
        bed = outDir.resolve(stem + "_targets.bed");
        System.out.println("Saving the target regions (BED) to a file: " + bed);
        Files.writeString(bed, plot.targetsBed());
    }

    /** The amplicons of target sequence {@code s} in one pool of the shared panel. */
    private static List<Pair> poolOf(List<Pair> pool, int s) {
        List<Pair> list = new ArrayList<>();
        for (Pair pcrx : pool) {
            if (pcrx.seq == s) {
                list.add(pcrx);
            }
        }
        return list;
    }

    /** Writes one pool to the panel report (sr1) and its primers to the primer list (srpanel). */
    private static void appendPanel(StringBuilder sr1, StringBuilder srpanel, String label, List<Pair> pcrpanel, String seq, RegionNames names) {
        if (pcrpanel.isEmpty()) {
            return;
        }
        sr1.append(label).append(":\n");
        sr1.append("Name\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)\n");
        for (Pair pcrx : pcrpanel) {
            srpanel.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\n");
            srpanel.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\n");
            sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
            sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.rCG)).append("\t").append(pcrx.rlc).append("\n");

            String amp = seq.substring(pcrx.fx5, 1 + pcrx.rx5);
            DNA2 dna2 = new DNA2(amp);
            double Tm = dna2.getTm65();
            double CG = dna2.getCG();
            sr1.append(names.design(pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp Tm=").append(String.format("%.1f", Tm)).append(" CG%=").append(String.format("%.1f", CG)).append(forcedMark(pcrx)).append("\n\n");
        }
    }

    /**
     * Heading of design region {@code h} in the primer list: its name, its
     * boundaries and, for a target merged into another region, which one.
     */
    private static String regionLine(RegionNames names, int[] exons, int[] cexons, int[] widened, int h) {
        int j = 2 * h;
        if (exons[j] < 0) {
            return names.target(h) + ": " + cexons[j] + "-" + cexons[j + 1] + " (" + (cexons[j + 1] - cexons[j] + 1) + "bp) join with region " + (-exons[j]);
        }
        return names.design(h) + ": " + exons[j] + "-" + exons[j + 1] + " " + (exons[j + 1] - exons[j] + 1) + "bp" + flankNote(widened, h);
    }

    /** Marks an amplicon line whose pair was forced into the panel. */
    private static String forcedMark(Pair pcrx) {
        return (pcrx.note == null) ? "" : " FORCED (" + pcrx.note + ")";
    }

    /** Notes a region rescued with widened primer-search flanks. */
    private static String flankNote(int[] widened, int h) {
        return (widened != null && h < widened.length && widened[h] > 0)
                ? " search flanks widened to " + widened[h] + " bp" : "";
    }

    private static int StrToInt(String str) {
        StringBuilder r = new StringBuilder();
        int z = 0;
        r.append(0);
        for (int i = 0; i < str.length(); i++) {
            char chr = str.charAt(i);
            if (chr > 47 && chr < 58) {
                r.append(chr);
                z++;
                if (z > 10) {
                    break;
                }
            }
            if (chr == '.' || chr == ',') {
                break;
            }
        }
        return (Integer.parseInt(r.toString()));
    }

    private static String PrimerReturn(String str) {
        String[] s1 = str.split("\t");
        String[] s2 = str.split(" ");
        int n1 = s1.length;
        int n2 = s2.length;
        String s = str;

        if (n1 > 1 && n2 < n1) {
            s = s1[1];
        }
        if (n2 > 1 && n1 < n2) {
            s = s2[1];
        }
        if (n1 == 1 && n2 == 1) {
            s = s1[0];
        }
        return dna.DNA(s);
    }

    public static String[] combineArrays(String[] array1, String[] array2) {
        String[] combinedArray = new String[array1.length + array2.length];
        System.arraycopy(array1, 0, combinedArray, 0, array1.length);
        System.arraycopy(array2, 0, combinedArray, array1.length, array2.length);
        return combinedArray;
    }

    static String stripQuotes(String s) {
        if (s == null) {
            return null;
        }
        s = s.trim();
        if ((s.startsWith("\"") && s.endsWith("\"")) || (s.startsWith("'") && s.endsWith("'"))) {
            if (s.length() >= 2) {
                return s.substring(1, s.length() - 1).trim();
            }
        }
        return s;
    }

    static String fileStem(Path p) {
        String name = p.getFileName().toString();
        int dot = name.lastIndexOf('.');
        return (dot > 0) ? name.substring(0, dot) : name;
    }

    /**
     * @param regions annotation of the exon pairs, for the panel map; null for FASTA
     */
    private static void RunSingleplexPanelDesign(LongIntHashMap map, String tagfile, String outpath, String seq, String name, String refsequence, String fastaseq, Boolean homology, String[] listprimers, int[] exons, int minpcr, int maxpcr, int minlen, int maxlen, int mintm, int maxtm, int minlc, int prlap, String ftail, String rtail, String e5, String e3, List<ExonPromoterExtractor.Region> regions) throws IOException {
        Path in = Paths.get(tagfile);
        String stem = fileStem(in);
        //String stem = fileStem(in).replaceFirst("\\.[0-9]+$", "");  // NG_029916.1

        Path outDir = (outpath == null || outpath.isBlank()) ? in.getParent() : Paths.get(outpath);
        if (outDir == null) {
            outDir = Paths.get(".");
        }
        // Path outDir = base.resolve(stem);

        Path primerlistfile = outDir.resolve(stem + "_primers.txt");
        Path pcrcolfile = outDir.resolve(stem + "_panels.txt");
        Path panelprimerlistfile = outDir.resolve(stem + "_panelprimers.txt");

        System.out.println(primerlistfile);
        System.out.println(pcrcolfile);
        System.out.println(panelprimerlistfile);

        int[] cexons = exons.clone();

        /*   ExonsDesigner pex = new ExonsDesigner(minpcr, maxpcr);
        List<ExonsDesigner.PcrFragment> fr = pex.buildPcrFragments(cexons);
        fr.forEach(System.out::println);
         */
        try {
            StringBuilder sr = new StringBuilder(100000);
            StringBuilder sr1 = new StringBuilder(100000);
            StringBuilder srpanel = new StringBuilder(100000);

            System.out.println("Running...");
            System.out.println("\nTarget file name: " + tagfile);

            sr.append("Target file name: ").append(tagfile).append("\n");
            sr.append("min PCR size=").append(minpcr).append("\n");
            sr.append("max PCR size=").append(maxpcr).append("\n");
            sr.append("min Tm=").append(mintm).append("\n");
            sr.append("max Tm=").append(maxtm).append("\n");
            sr.append("min length=").append(minlen).append("\n");

            if (e3.length() > 0) {
                sr.append("3-end=").append(e3).append("\n");
            }
            if (e5.length() > 0) {
                sr.append("5-end=").append(e5).append("\n");
            }
            if (ftail.length() > 0) {
                sr.append("Forward tail=").append(ftail).append("\n");
            }
            if (rtail.length() > 0) {
                sr.append("Reverse tail=").append(rtail).append("\n");
            }

            PrimerDesign pd = new PrimerDesign(map, seq, name, minlen, maxlen, mintm, maxtm, minlc, exons, listprimers, refsequence, fastaseq, homology);

            List<Pair> pcrcol = pd.runDesign(ftail, rtail, e5, e3, prlap, minpcr, maxpcr);

            RegionNames names = new RegionNames(exons, regions);
            PanelsCollector panels = new PanelsCollector(pcrcol, pd.getAllPairs());
            panels.setRegionNames(List.of(names));
            panels.CombinePairsPanel(cexons, 1);
            List<Pair> pcrpanel = panels.getPanel1();

            if (!pcrpanel.isEmpty()) {
                sr1.append("Panel:\n");
                sr1.append("Name\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)\n");
                for (Pair pcrx : pcrpanel) {
                    srpanel.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\n");
                    srpanel.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\n");
                    sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
                    sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.rCG)).append("\t").append(pcrx.rlc).append("\n");

                    String amp = seq.substring(pcrx.fx5, 1 + pcrx.rx5);
                    DNA2 dna2 = new DNA2(amp);
                    double Tm = dna2.getTm65();
                    double CG = dna2.getCG();
                    sr1.append(names.design(pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp Tm=").append(String.format("%.1f", Tm)).append(" CG%=").append(String.format("%.1f", CG)).append(forcedMark(pcrx)).append("\n");
                    // sr1.append(">\n").append(amp).append("\n\n");
                    sr1.append(amp).append("\n\n");
                }
            }

            String notes = panels.renderNotes(0);
            if (!notes.isEmpty()) {
                sr1.append(notes);
                System.out.println(notes);
            }

            if (!pcrpanel.isEmpty()) {
                CoverageReport cov = new CoverageReport(cexons, regions);
                List<int[]> iv = new ArrayList<>();
                for (Pair pcrx : pcrpanel) {
                    iv.add(CoverageReport.interval(pcrx.fx5, 1 + pcrx.rx5));
                }
                cov.addPool("Panel", iv);
                String covtxt = cov.render();
                if (!covtxt.isEmpty()) {
                    sr1.append(covtxt);
                    System.out.println(covtxt);
                }
            }

            PrimersCollector[] fPrimersList = pd.getForwardPrimers();
            PrimersCollector[] rPrimersList = pd.getReversePrimers();
            int[] widened = pd.getWidenedFlanks();
            int h = -1;
            sr.append("\nName\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)");
            for (int j = 0; j < exons.length - 1; j += 2) {
                h++;
                String region = regionLine(names, exons, cexons, widened, h);
                sr.append("\n").append(region).append("\n");
                System.out.println("\n" + region);

                PrimersCollector fPrimersList1 = fPrimersList[h];
                if (fPrimersList1.Amount() > 0) {
                    StringBuilder sr2 = new StringBuilder(100000);
                    double[] Tm = fPrimersList1.getTms();
                    double[] CG = fPrimersList1.getCGs();
                    int[] ln = fPrimersList1.getPrimerLengths();          // length
                    int[] lc = fPrimersList1.getPrimerLC();               // Linguistic_Complexity
                    // int[] x1 = fPrimersList1.getPrimerLocations();        // location x1
                    String[] primer = fPrimersList1.getPrimer();          // primer sequence
                    String[] primername = fPrimersList1.getpPrimerName(); // name

                    for (int i = 0; i < fPrimersList1.Amount(); i++) {
                        sr2.append(primername[i]).append("\t").append(primer[i]).append("\t").append(ln[i]).append("\t").append(String.format("%.1f", Tm[i])).append("\t").append(String.format("%.1f", CG[i])).append("\t").append(lc[i]).append("\n");
                    }
                    System.out.println(sr2);
                    sr.append(sr2);
                }

                PrimersCollector rPrimersList1 = rPrimersList[h];
                if (rPrimersList1.Amount() > 0) {
                    StringBuilder sr2 = new StringBuilder(100000);
                    double[] Tm = rPrimersList1.getTms();
                    double[] CG = rPrimersList1.getCGs();
                    int[] ln = rPrimersList1.getPrimerLengths();          // length
                    int[] lc = rPrimersList1.getPrimerLC();               // Linguistic_Complexity
                    //  int[] x1 = rPrimersList1.getPrimerLocations();        // location x1
                    String[] primer = rPrimersList1.getPrimer();          // primer sequence
                    String[] primername = rPrimersList1.getpPrimerName(); //name
                    for (int i = 0; i < rPrimersList1.Amount(); i++) {
                        sr2.append(primername[i]).append("\t").append(primer[i]).append("\t").append(ln[i]).append("\t").append(String.format("%.1f", Tm[i])).append("\t").append(String.format("%.1f", CG[i])).append("\t").append(lc[i]).append("\n");
                    }
                    System.out.println(sr2);
                    sr.append(sr2);
                }
            }
/*
            try (FileWriter fileWriter = new FileWriter(primerlistfile.toFile())) {
                System.out.println("Saving the primer list report to a file: " + primerlistfile);
                fileWriter.write(sr.toString());
            }

            if (!pcrpanel.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(panelprimerlistfile.toFile())) {
                    System.out.println("Saving panels PCR primers combinations report to a file: " + panelprimerlistfile);
                    fileWriter.write(srpanel.toString());
                }
            }
*/
            PanelPlot plot = new PanelPlot(name, seq.length(), cexons, exons, regions, pd.getWidenedFlanks(), List.of(pcrpanel), false);
            if (!pcrcol.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(pcrcolfile.toFile())) {
                    System.out.println("Saving PCR primers combinations report to a file: " + pcrcolfile);
                    fileWriter.write(sr1.toString());
                }
                writePlot(outDir, stem, plot);
            }
            writeUncovered(outDir, stem, plot.uncoveredFasta(seq, maxpcr));

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    static LongIntHashMap GenomeReference(List<String> genomefiles) throws IOException {
        List<String> seqs = new ArrayList<>();  // Initialize the list

        for (String files : genomefiles) {
            List<FastaIO.FastaRecord> records = FastaIO.readAllToList(java.util.Collections.singletonList(files));
            if (records.isEmpty()) {
                continue;
            }
            for (FastaIO.FastaRecord rec : records) {
                if (rec == null || rec.sequence == null || rec.sequence.isBlank()) {
                    continue;
                }
                seqs.add(rec.sequence);
            }
        }

        MaskingSequences ms = new MaskingSequences();
        return ms.mask(seqs.toArray(String[]::new), 18);//int GENOME_KMER_SIZE = 18;
    }

    private static Path validatePath(String line, String key, boolean isDirectory, boolean allowMissing) {
        if (!line.toLowerCase().contains(key.toLowerCase())) {
            return null;
        }
        int i = line.toLowerCase().indexOf(key.toLowerCase());
        String pathStr = line.substring(i + key.length()).trim();
        pathStr = stripQuotes(pathStr);
        if (pathStr.isEmpty()) {
            return null;
        }
        try {
            Path path = Paths.get(pathStr);
            if (!Files.exists(path)) {
                if (allowMissing) {
                    return path;
                }
                System.err.println("Error: Path does not exist: " + path);
                return null;
            }
            if (isDirectory && !Files.isDirectory(path)) {
                System.err.println("Error: Path is not a directory: " + path);
                return null;
            }
            if (!isDirectory && !Files.isRegularFile(path)) {
                System.err.println("Error: Path is not a file: " + path);
                return null;
            }
            if (!Files.isReadable(path)) {
                System.err.println("Error: No read permission: " + path);
                return null;
            }
            return path;
        } catch (Exception e) {
            System.err.println("Error processing path: " + e.getMessage());
            return null;
        }
    }
}
