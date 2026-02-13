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
import java.util.Set;

public class apanel {

    // ========================== CONFIG CLASS ==========================

    /**
     * Holds all primer-design parameters parsed from the configuration file.
     * Replaces 15+ loose local variables with a single, passable object.
     */
    static class PrimerConfig {
        int minPcr = 60;
        int maxPcr = 600;
        int minTm  = 60;
        int maxTm  = 62;
        int minLen = 18;
        int maxLen = 25;
        int prLap  = 12;   // 0-18
        int minLc  = 70;
        String e5     = "n";
        String e3     = "w";
        String fTail  = "";
        String rTail  = "";
        String outPath = null;

        boolean homology  = false;
        boolean multiplex = true;
    }

    // ========================== HELPERS ==========================

    /** Clamp value into [min, max]. */
    private static int clamp(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }

    /**
     * Extract the value portion after the FIRST '=' on the line.
     * Returns "" if no '=' is found.
     */
    private static String valueAfterEquals(String line) {
        int idx = line.indexOf('=');
        return (idx >= 0) ? line.substring(idx + 1).trim() : "";
    }

    /**
     * Parse an integer parameter from a config line that starts with the given key.
     * Returns null if the line doesn't contain the key.
     */
    private static Integer parseIntParam(String line, String key) {
        if (!line.contains(key)) {
            return null;
        }
        System.out.println(line);
        return StrToInt(valueAfterEquals(line));
    }

    /**
     * Build the concatenation of all sequences EXCEPT the one at index {@code excludeIdx}.
     */
    private static String buildOtherSequences(List<String> sequences, int excludeIdx) {
        int totalLen = 0;
        for (int j = 0; j < sequences.size(); j++) {
            if (j == excludeIdx) continue;
            String s = sequences.get(j);
            if (s != null) totalLen += s.length();
        }
        StringBuilder sb = new StringBuilder(totalLen);
        for (int j = 0; j < sequences.size(); j++) {
            if (j == excludeIdx) continue;
            String s = sequences.get(j);
            if (s != null) sb.append(s);
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
     * Read and concatenate reference FASTA file(s) into a single cleaned sequence string.
     */
    private static String readReferenceSequences(List<String> refFiles) {
        if (refFiles.isEmpty()) return "";

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
     * Load sequence + name + exons from a tagfile, with FASTA fallback.
     * Returns null if no valid sequence is found.
     */
    private static class SeqRecord {
        final String sequence;
        final String name;
        final int[] exons;
        SeqRecord(String sequence, String name, int[] exons) {
            this.sequence = sequence;
            this.name = name;
            this.exons = exons;
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
            return new SeqRecord(sequence, name, exons);
        }

        // FASTA fallback — use last valid record
        List<FastaIO.FastaRecord> records =
                FastaIO.readAllToList(java.util.Collections.singletonList(tagFile));
        for (FastaIO.FastaRecord rec : records) {
            if (rec != null && rec.sequence != null && !rec.sequence.isBlank()) {
                String recName = (rec.id != null && !rec.id.isBlank()) ? rec.id : name;
                int[] exons = new int[]{0, rec.sequence.length()};
                // Note: returns LAST valid record (preserving original behavior)
                sequence = rec.sequence;
                name = recName;
                return new SeqRecord(sequence, name, exons);
            }
        }
        return null; // nothing found
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

        try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String cline = line;       // preserve original case for paths
                line = line.toLowerCase();  // lowercase for key matching

                // --- Path-based parameters ---
                Path folderPath = validatePath(cline, "folder_path=", true, false);
                if (folderPath != null) {
                    String s = cline.substring(cline.indexOf("folder_path=") + "folder_path=".length()).trim();
                    System.out.println("Target_Folder_path=" + s);
                    String[] files = new BatchReadFiles().BatchReadFiles(new String[]{s, "*.*"});
                    tagfiles.addAll(List.of(files));
                }

                Path genomePath = validatePath(cline, "genome_path=", true, false);
                if (genomePath != null) {
                    String s = cline.substring(cline.indexOf("genome_path=") + "genome_path=".length()).trim();
                    System.out.println("Genome_path=" + genomePath);
                    String[] files = new BatchReadFiles().BatchReadFiles(new String[]{s, "*.*"});
                    genomefiles.addAll(List.of(files));
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
                    System.out.println("Reference_file_path= " + refFile);  // BUG FIX: was printing outDir
                    reffiles.add(refFile.toString());
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

                // --- String parameters (use startsWith to avoid partial matches) ---
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
                if ((v = parseIntParam(line, "minpcr=")) != null) cfg.minPcr = clamp(v, 30, 5000);
                if ((v = parseIntParam(line, "maxpcr=")) != null) cfg.maxPcr = clamp(v, 30, 50000);
                if ((v = parseIntParam(line, "minlen=")) != null) cfg.minLen = clamp(v, 12, 80);
                if ((v = parseIntParam(line, "maxlen=")) != null) cfg.maxLen = clamp(v, 12, 100);
                if ((v = parseIntParam(line, "mintm="))  != null) cfg.minTm  = clamp(v, 40, 75);
                if ((v = parseIntParam(line, "maxtm="))  != null) cfg.maxTm  = clamp(v, 40, 80);
            }
        } catch (IOException e) {
            System.err.println("Error reading config file '" + infile + "': " + e.getMessage());
        }
        return cfg;
    }

    // ========================== SINGLE-PLEX ==========================

    private static void runSingleplexAll(Set<Long> map, List<String> tagfiles,
                                         PrimerConfig cfg, String refsequence,
                                         String[] primers) throws Exception {
        final String ref = "";
        for (String tagfile : tagfiles) {
            ExonPromoterExtractor ex = new ExonPromoterExtractor(tagfile);
            String sequence = ex.getSequence();
            String name = ex.getRecordId();

            // FASTA fallback
            if (sequence == null || sequence.isBlank()) {
                List<FastaIO.FastaRecord> records =
                        FastaIO.readAllToList(java.util.Collections.singletonList(tagfile));
                if (records.isEmpty()) continue;

                for (FastaIO.FastaRecord rec : records) {
                    if (rec == null || rec.sequence == null || rec.sequence.isBlank()) continue;
                    String recName = (rec.id != null && !rec.id.isBlank()) ? rec.id : name;
                    int[] exons = new int[]{0, rec.sequence.length()};
                    RunSingleplexPanelDesign(map, tagfile, cfg.outPath, rec.sequence, recName,
                            refsequence, ref, cfg.homology, primers, exons,
                            cfg.minPcr, cfg.maxPcr, cfg.minLen, cfg.maxLen,
                            cfg.minTm, cfg.maxTm, cfg.minLc, cfg.prLap,
                            cfg.fTail, cfg.rTail, cfg.e5, cfg.e3);
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
                    cfg.fTail, cfg.rTail, cfg.e5, cfg.e3);
        }
    }

    // ========================== MULTIPLEX ==========================

    private static void runMultiplexAll(Set<Long> map, List<String> tagfiles,
                                        PrimerConfig cfg, String refsequence,
                                        String[] primers) throws Exception {
        List<String> sequences     = new ArrayList<>();
        List<String> nameSequences = new ArrayList<>();
        List<int[]>  exonsSequences = new ArrayList<>();

        // Phase 1: load all sequences
        for (String tagfile : tagfiles) {
            SeqRecord rec = loadSequenceRecord(tagfile);
            if (rec == null) continue;
            sequences.add(rec.sequence);
            nameSequences.add(rec.name);
            exonsSequences.add(rec.exons);
        }

        // Phase 2: design primers file-by-file
        int nseq = sequences.size();
        for (int i = 0; i < nseq; i++) {
            String sequence = sequences.get(i);
            if (sequence == null || sequence.isEmpty()) continue;

            String name    = nameSequences.get(i);
            int[]  exons   = exonsSequences.get(i);
            String tagfile = tagfiles.get(i);
            String otherSeqs = buildOtherSequences(sequences, i);

            String[] newPrimers = RunMultiplexPanelDesign(map, tagfile, cfg.outPath,
                    sequence, name, refsequence, otherSeqs,
                    cfg.homology, primers, exons,
                    cfg.minPcr, cfg.maxPcr, cfg.minLen, cfg.maxLen,
                    cfg.minTm, cfg.maxTm, cfg.minLc, cfg.prLap,
                    cfg.fTail, cfg.rTail, cfg.e5, cfg.e3);
            primers = combineArrays(primers, newPrimers);
        }
    }

    // ========================== MAIN ==========================

    public static void main(String[] args) throws IOException, Exception {
        if (args.length == 0) {
            System.err.println("Usage: <program> <config-file>");
            return;
        }

        String infile = args[0];
        System.out.println("Current Directory: " + System.getProperty("user.dir"));
        System.out.println("Command-line arguments: " + infile);

        // Mutable lists populated by parseConfig
        List<String> tagfiles    = new ArrayList<>();
        List<String> reffiles    = new ArrayList<>();
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
        Set<Long> map = GenomeReference(genomefiles);
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

    private static String[] RunMultiplexPanelDesign(Set<Long> map, String tagfile, String outpath, String seq, String name, String refsequence, String fastaseq, Boolean homology, String[] listprimers, int[] exons, int minpcr, int maxpcr, int minlen, int maxlen, int mintm, int maxtm, int minlc, int prlap, String ftail, String rtail, String e5, String e3) throws IOException {
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

        String[] listprimers2 = null;
        int[] cexons = exons.clone();

        /*   ExonsDesigner pex = new ExonsDesigner(minpcr, maxpcr);
        List<ExonsDesigner.PcrFragment> fr = pex.buildPcrFragments(cexons);
        fr.forEach(System.out::println);
         */
        try {
            StringBuilder sr = new StringBuilder(100000);
            StringBuilder sr1 = new StringBuilder(100000);
            StringBuilder srpanel = new StringBuilder(100000);

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
            List<Pair> pcrcol = pd.RunDesign(ftail, rtail, e5, e3, prlap, minpcr, maxpcr);

            PanelsCollector panels = new PanelsCollector(pcrcol);
            panels.CombinePairsPanel(exons.length);
            List<Pair> pcrpanel1 = panels.getPanel1();
            List<Pair> pcrpanel2 = panels.getPanel2();
            ArrayList<String> stringList = new ArrayList<>();

            if (!pcrpanel1.isEmpty()) {
                sr1.append("Panel1:\n");
                sr1.append("Name\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)\n");
                for (Pair pcrx : pcrpanel1) {
                    stringList.add(pcrx.fprimer);
                    stringList.add(pcrx.rprimer);
                    srpanel.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\n");
                    srpanel.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\n");
                    sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
                    sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.rCG)).append("\t").append(pcrx.rlc).append("\n");

                    String amp = seq.substring(pcrx.fx5, pcrx.rx5);
                    DNA2 dna2 = new DNA2(amp);
                    double Tm = dna2.getTm65();
                    double CG = dna2.getCG();
                    sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp Tm=").append(String.format("%.1f", Tm)).append(" CG%=").append(String.format("%.1f", CG)).append("\n\n");
                }
            }
            if (!pcrpanel2.isEmpty()) {
                sr1.append("Panel2:\n");
                sr1.append("Name\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)\n");
                for (Pair pcrx : pcrpanel2) {
                    stringList.add(pcrx.fprimer);
                    stringList.add(pcrx.rprimer);
                    srpanel.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\n");
                    srpanel.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\n");
                    sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
                    sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.rCG)).append("\t").append(pcrx.rlc).append("\n");
                    //   sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp\n\n");

                    String amp = seq.substring(pcrx.fx5, pcrx.rx5);
                    DNA2 dna2 = new DNA2(amp);
                    double Tm = dna2.getTm65();
                    double CG = dna2.getCG();
                    sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp Tm=").append(String.format("%.1f", Tm)).append(" CG%=").append(String.format("%.1f", CG)).append("\n\n");

                }
            }
            listprimers2 = stringList.toArray(String[]::new);
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
            PrimersCollector[] fPrimersList = pd.getpForwardPrimers();
            PrimersCollector[] rPrimersList = pd.getpReversePrimers();
            int h = -1;
            sr.append("\nName\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)");
            for (int j = 0; j < exons.length - 1; j += 2) {
                h++;
                if (exons[j] < 0) {
                    sr.append("\n").append("exon:").append(h + 1).append(" ").append(cexons[j]).append("-").append(cexons[j + 1]).append(" (").append(cexons[j + 1] - cexons[j] + 1).append("bp) join with ").append(-exons[j]).append("\n");
                    System.out.println("\nexon:" + (h + 1) + " " + cexons[j] + "-" + cexons[j + 1] + " (" + (cexons[j + 1] - cexons[j] + 1) + "bp) join with " + (-exons[j]));
                } else {
                    sr.append("\n").append("exon:").append(h + 1).append(" ").append(exons[j]).append("-").append(exons[j + 1]).append(" ").append(exons[j + 1] - exons[j] + 1).append("bp\n");
                    System.out.println("\nexon:" + (h + 1) + " " + exons[j] + "-" + exons[j + 1] + " (" + (exons[j + 1] - exons[j] + 1) + "bp)");
                }

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

            if (!pcrpanel1.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(panelprimerlistfile.toFile())) {
                    System.out.println("Saving panels PCR primers combinations report to a file: " + panelprimerlistfile);
                    fileWriter.write(srpanel.toString());
                }
            }

            if (!pcrcol.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(pcrcolfile.toFile())) {
                    System.out.println("Saving PCR primers combinations report to a file: " + pcrcolfile);
                    fileWriter.write(sr1.toString());
                }
            }
            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds\n\n");

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
        return listprimers2;

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

    private static void RunSingleplexPanelDesign(Set<Long> map, String tagfile, String outpath, String seq, String name, String refsequence, String fastaseq, Boolean homology, String[] listprimers, int[] exons, int minpcr, int maxpcr, int minlen, int maxlen, int mintm, int maxtm, int minlc, int prlap, String ftail, String rtail, String e5, String e3) throws IOException {
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

            List<Pair> pcrcol = pd.RunDesign(ftail, rtail, e5, e3, prlap, minpcr, maxpcr);

            PanelsCollector panels = new PanelsCollector(pcrcol);
            panels.CombinePairsPanel(exons.length);
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
                    sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp Tm=").append(String.format("%.1f", Tm)).append(" CG%=").append(String.format("%.1f", CG)).append("\n");
                    // sr1.append(">\n").append(amp).append("\n\n");
                    sr1.append(amp).append("\n\n");
                }
            }

            PrimersCollector[] fPrimersList = pd.getpForwardPrimers();
            PrimersCollector[] rPrimersList = pd.getpReversePrimers();
            int h = -1;
            sr.append("\nName\tSequence\tLength\tTm(°C)\tGC(%)\tLinguistic_Complexity(%)");
            for (int j = 0; j < exons.length - 1; j += 2) {
                h++;
                if (exons[j] < 0) {
                    sr.append("\n").append("exon:").append(h + 1).append(" ").append(cexons[j]).append("-").append(cexons[j + 1]).append(" (").append(cexons[j + 1] - cexons[j] + 1).append("bp) join with ").append(-exons[j]).append("\n");
                    System.out.println("\nexon:" + (h + 1) + " " + cexons[j] + "-" + cexons[j + 1] + " (" + (cexons[j + 1] - cexons[j] + 1) + "bp) join with " + (-exons[j]));
                } else {
                    sr.append("\n").append("exon:").append(h + 1).append(" ").append(exons[j]).append("-").append(exons[j + 1]).append(" ").append(exons[j + 1] - exons[j] + 1).append("bp\n");
                    System.out.println("\nexon:" + (h + 1) + " " + exons[j] + "-" + exons[j + 1] + " (" + (exons[j + 1] - exons[j] + 1) + "bp)");
                }

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
            if (!pcrcol.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(pcrcolfile.toFile())) {
                    System.out.println("Saving PCR primers combinations report to a file: " + pcrcolfile);
                    fileWriter.write(sr1.toString());
                }
            }

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    static Set<Long> GenomeReference(List<String> genomefiles) throws IOException {
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
        return ms.mask(seqs.toArray(String[]::new), 19);
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
