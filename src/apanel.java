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

    public static void main(String[] args) throws IOException, Exception {
        if (args.length > 0) {
            String[] primers = new String[0];
            String primerfile = "";
            String infile = args[0];
            for (String arg : args) {
                if (args.length > 0) {
                    infile = arg;
                    break;
                }
            }

            List<String> Sequences = new ArrayList<>();
            List<String> NameSequences = new ArrayList<>();
            List<int[]> ExonsSequences = new ArrayList<>();
            List<String> tagfiles = new ArrayList<>();
            List<String> tagprimers = new ArrayList<>();
            List<String> reffiles = new ArrayList<>();
            List<String> genomefiles = new ArrayList<>();

            int minpcr = 60;
            int maxpcr = 600;
            int mintm = 60;
            int maxtm = 62;
            int minlen = 18;
            int maxlen = 25;
            int prlap = 12; //0-18
            int minlc = 70;
            String e5 = "n";
            String e3 = "w";
            String ftail = "";
            String rtail = "";
            String outpath = null;

            Boolean homology = false;
            Boolean multiplex = true;

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
                String line;

                while ((line = br.readLine()) != null) {
                    String cline = line;
                    line = line.toLowerCase();

//Path file = validatePath(line, "input_file=", false);  // file
//Path dir = validatePath(line, "reference_path=", true); // directory
                    Path outDir = validatePath(cline, "folder_path=", true);
                    if (outDir != null) {
                        String s = cline.substring(12);
                        System.out.println("Folder_path=" + s);    // Reading target file(s), one by one
                        String[] files = new BatchReadFiles().BatchReadFiles(new String[]{s, "*.*"});
                        tagfiles.addAll(List.of(files));
                    }

                    Path inputFile = validatePath(cline, "target_path=", false);
                    if (inputFile != null) {
                        System.out.println("Input file= " + inputFile);
                    //    tagfiles.addAll(List.of(inputFile.toString()));
                        tagfiles.add(inputFile.toString()); 
                    }

                    outDir = validatePath(cline, "genome_path=", true);
                    if (outDir != null) {
                        System.out.println("Genome_path=\"" + outDir);
                        String s = cline.substring(12);
                        String[] files = new BatchReadFiles().BatchReadFiles(new String[]{s, "*.*"});
                        genomefiles = List.of(files);
                    }

                    outDir = validatePath(cline, "folder_out=", true);
                    if (outDir != null) {
                        int i = cline.indexOf("folder_out=");
                        outpath = cline.substring(i + "folder_out=".length()).trim();
                        outpath = stripQuotes(outpath);
                        outDir = OutDirUtil.prepareOutputDir(outpath);
                        System.out.println("Output dir ready: " + outDir);
                        outpath = outDir.toString();
                    }

                    outDir = validatePath(cline, "reference_path=", true);
                    if (outDir != null) {
                        System.out.println("Reference_path= " + outDir);
                        reffiles.add(outDir.toString());
                    }

                    if (line.contains("homology=true")) {
                        homology = true;
                        System.out.println("Designing common primers only based on shared sequences between different files.");
                    }

                    if (line.contains("multiplex=false")) {
                        multiplex = false;
                        System.out.println("");
                    }
                    if (line.contains("3end=")) {
                        e3 = line.substring(5);
                        System.out.println("3End=" + e3);
                    }
                    if (line.contains("5end=")) {
                        e5 = line.substring(5);
                        System.out.println("5End=" + e5);
                    }
                    if (line.contains("target_primers=")) {
                        String s = line.substring(15);
                        primerfile = s;
                        System.out.println("Target_primers=" + s);
                    }
                    if (line.contains("forwardtail=")) {
                        String s = dna.DNA(line.substring(12));
                        ftail = s;
                        System.out.println("ForwardTail=" + s);
                    }
                    if (line.contains("reversetail=")) {
                        String s = dna.DNA(line.substring(12));
                        rtail = s;
                        System.out.println("ReverseTail=" + s);
                    }
                    if (line.contains("minpcr=")) {
                        int h = StrToInt(line.substring(7));
                        System.out.println(line);
                        if (h < 30) {
                            h = 30;
                        }
                        if (h > 5000) {
                            h = 5000;
                        }
                        minpcr = h;
                    }
                    if (line.contains("maxpcr=")) {
                        int h = StrToInt(line.substring(7));
                        System.out.println(line);
                        if (h < 30) {
                            h = 30;
                        }
                        if (h > 50000) {
                            h = 50000;
                        }
                        maxpcr = h;
                    }
                    if (line.contains("minlen=")) {
                        int h = StrToInt(line.substring(6));
                        System.out.println(line);
                        if (h < 12) {
                            h = 12;
                        }
                        if (h > 80) {
                            h = 80;
                        }
                        minlen = h;
                    }
                    if (line.contains("maxlen=")) {
                        int h = StrToInt(line.substring(6));
                        System.out.println(line);
                        if (h < 12) {
                            h = 12;
                        }
                        if (h > 100) {
                            h = 100;
                        }
                        maxlen = h;
                    }
                    if (line.contains("mintm=")) {
                        int h = StrToInt(line.substring(5));
                        System.out.println(line);
                        if (h < 40) {
                            h = 40;
                        }
                        if (h > 75) {
                            h = 75;
                        }
                        mintm = h;
                    }
                    if (line.contains("maxtm=")) {
                        int h = StrToInt(line.substring(5));
                        System.out.println(line);
                        if (h < 40) {
                            h = 40;
                        }
                        if (h > 80) {
                            h = 80;
                        }
                        maxtm = h;
                    }
                }
            } catch (IOException e) {
            }

// Reading file with primers list
            if (primerfile.length() > 1) {
                try (BufferedReader bufferedReader = new BufferedReader(new FileReader(primerfile))) {
                    String line;
                    int k = 0;
                    int n = -1;
                    while ((line = bufferedReader.readLine()) != null) {
                        line = line.trim().toLowerCase();
                        if (!line.isBlank()) {
                            k++;
                            tagprimers.add(line);
                        }
                    }
                    primers = new String[k];
                    for (String p : tagprimers) {
                        p = PrimerReturn(p);
                        if (p.length() > 8) {
                            n++;
                            primers[n] = p;
                        }
                    }
                } catch (IOException e) {
                }
            }

// Reading reference FASTA file(s), one by one
            String refsequence = ""; // FASTA sequence for repeated blocks identification
            if (!reffiles.isEmpty()) {
                for (String reffile : reffiles) {
                    try (BufferedReader bufferedReader = new BufferedReader(
                            new InputStreamReader(
                                    reffile.endsWith(".gz")
                                    ? new java.util.zip.GZIPInputStream(new FileInputStream(reffile))
                                    : new FileInputStream(reffile),
                                    java.nio.charset.StandardCharsets.UTF_8), 1 << 16)) {

                        String line;
                        StringBuilder sequenceBuilder = new StringBuilder(1 << 20);

                        while ((line = bufferedReader.readLine()) != null) {
                            if (line.isEmpty()) {
                                continue;
                            }
                            if (line.charAt(0) == '>') {
                            } else {
                                sequenceBuilder.append(line);
                            }
                        }
                        refsequence = refsequence + dna.DNA(sequenceBuilder.toString());
                    }
                }
            }

            int nseq = 0;

//*********reading genome files for hash-map
            Set<Long> map = GenomeReference(genomefiles);

//*********single-plex************* special task for exons primer 
            if (!multiplex) {
                final String ref = ""; // keep if the design API expects it

                for (String tagfile : tagfiles) {
                    ExonPromoterExtractor ex = new ExonPromoterExtractor(tagfile);
                    String sequence = ex.getSequence();
                    String name = ex.getRecordId(); // or ex.getAccession()

                    // ---- FASTA fallback if extractor didn't yield sequence ----
                    if (sequence == null || sequence.isBlank()) {
                        List<FastaIO.FastaRecord> records = FastaIO.readAllToList(java.util.Collections.singletonList(tagfile));
                        if (records.isEmpty()) {
                            continue;
                        }
                        for (FastaIO.FastaRecord rec : records) {
                            if (rec == null || rec.sequence == null || rec.sequence.isBlank()) {
                                continue;
                            }
                            name = (rec.id != null && !rec.id.isBlank()) ? rec.id : name;
                            int[] exons = new int[]{0, rec.sequence.length()};
                            RunSingleplexPanelDesign(map, tagfile, outpath, rec.sequence, name, refsequence, ref, homology, primers, exons, minpcr, maxpcr, minlen, maxlen, mintm, maxtm, minlc, prlap, ftail, rtail, e5, e3);
                        }
                        // IMPORTANT: avoid running the non-FASTA path after processing FASTA
                        continue;
                    }

                    // ---- GenBank / extracted regions path ----
                    List<ExonPromoterExtractor.Region> regs = ex.getRegions();
                    int[] exons = toExonArray(regs, sequence.length());
                    RunSingleplexPanelDesign(map, tagfile, outpath, sequence, name, refsequence, ref, homology, primers, exons, minpcr, maxpcr, minlen, maxlen, mintm, maxtm, minlc, prlap, ftail, rtail, e5, e3);
                }
                return;
            }

//********************** multiplex PCR - panel 
            for (String tagfile : tagfiles) {
                nseq++;
                ExonPromoterExtractor ex1 = new ExonPromoterExtractor(tagfile);
                List<ExonPromoterExtractor.Region> regs = ex1.getRegions();
                String sequence = ex1.getSequence();
                String name = ex1.getRecordId();
                int[] exons = regs.stream().flatMapToInt(r -> java.util.stream.IntStream.of(r.start, r.end)).toArray();

                if (sequence == null || sequence.isBlank()) {
                    List<FastaIO.FastaRecord> records = FastaIO.readAllToList(java.util.Collections.singletonList(tagfile));
                    if (records.isEmpty()) {
                        continue;
                    }
                    for (FastaIO.FastaRecord rec : records) {
                        if (rec == null || rec.sequence == null || rec.sequence.isBlank()) {
                            continue;
                        }
                        name = (rec.id != null && !rec.id.isBlank()) ? rec.id : name;
                        exons = new int[]{0, rec.sequence.length()};
                        sequence = rec.sequence;
                    }
                }
                ExonsSequences.add(exons);
                Sequences.add(sequence);
                NameSequences.add(name);
            }

//Running file by file            
            for (int i = 0; i < nseq; i++) {
                int[] exons = ExonsSequences.get(i);
                String sequence = Sequences.get(i);
                String name = NameSequences.get(i);
                String tagfile = tagfiles.get(i);

                if (!sequence.isEmpty()) {

                    String ref = "";
                    int totalLen = ref.length();
                    for (int j = 0; j < nseq; j++) {
                        if (j == i) {
                            continue;
                        }
                        String s = Sequences.get(j);
                        if (s != null) {
                            totalLen += s.length();
                        }
                    }
                    StringBuilder sb = new StringBuilder(totalLen);
                    sb.append(ref);
                    for (int j = 0; j < nseq; j++) {
                        if (j == i) {
                            continue;
                        }
                        String s = Sequences.get(j);
                        if (s != null) {
                            sb.append(s);
                        }
                    }
                    String[] newlistprimers = RunMultiplexPanelDesign(map, tagfile, outpath, sequence, name, refsequence, sb.toString(), homology, primers, exons, minpcr, maxpcr, minlen, maxlen, mintm, maxtm, minlc, prlap, ftail, rtail, e5, e3);
                    // Combine the primers with new set of primers
                    primers = combineArrays(primers, newlistprimers);
                }

            }

        }
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
        String stem = fileStem(in).replaceFirst("\\.[0-9]+$", "");  // NG_029916.1

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

            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds\n\n");
            try (FileWriter fileWriter = new FileWriter(primerlistfile.toFile())) {
                System.out.println("Saving the primer list report to a file: " + primerlistfile);
                fileWriter.write(sr.toString());
                //    fileWriter.write("\nTime taken: " + duration + " seconds\n\n");
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
        String stem = fileStem(in).replaceFirst("\\.[0-9]+$", "");  // NG_029916.1

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

            long startTime = System.nanoTime();
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

            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds\n\n");

            try (FileWriter fileWriter = new FileWriter(primerlistfile.toFile(), true)) {
                System.out.println("Saving the primer list report to a file: " + primerlistfile);
                fileWriter.write(sr.toString());
                //    fileWriter.write("\nTime taken: " + duration + " seconds\n\n");
            }

            if (!pcrpanel.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(panelprimerlistfile.toFile(), true)) {
                    System.out.println("Saving panels PCR primers combinations report to a file: " + panelprimerlistfile);
                    fileWriter.write(srpanel.toString());
                }
            }

            if (!pcrcol.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(pcrcolfile.toFile(), true)) {
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
        return ms.Mask(seqs.toArray(String[]::new), 21);
    }

    private static Path validatePath(String line, String key, boolean isDirectory) {
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
