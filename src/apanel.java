import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class apanel {

    public static void main(String[] args) {
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

            List<String> tagfiles = new ArrayList<>();
            List<String> tagprimers = new ArrayList<>();

            int minpcr = 60;
            int maxpcr = 600;
            int mintm = 60;
            int maxtm = 62;
            int minlen = 18;
            int maxlen = 25;
            int prlap = 12; //12-18
            int minlc = 78;
            String e5 = "n";
            String e3 = "w";
            String ftail = "";
            String rtail = "";

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    line = line.toLowerCase();

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
                    if (line.contains("target_path=")) {
                        String s = line.substring(12);
                        tagfiles.add(s);
                        System.out.println("Target_path=" + s);
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
                    if (line.contains("minlc=")) {
                        int h = StrToInt(line.substring(3));
                        System.out.println(line);
                        if (h < 20) {
                            h = 20;
                        }
                        if (h > 90) {
                            h = 90;
                        }
                        minlc = h;
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
                        if (h > 5000) {
                            h = 5000;
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

            // Reading target file(s), one by one
            for (String tagfile : tagfiles) {
                try (BufferedReader bufferedReader = new BufferedReader(new FileReader(tagfile))) {
                    String line;
                    String name = "";

                    StringBuilder exonsBuilder = new StringBuilder();
                    StringBuilder sequenceBuilder = new StringBuilder();

                    boolean isOriginSection = false;
                    int readingExons = 0;

                    while ((line = bufferedReader.readLine()) != null) {
                        line = line.trim().toLowerCase();
                        if (line.contains("locus")) {
                            String[] parts = line.split(" ");
                            int n = 0;
                            for (String part : parts) {
                                String a = part.trim();
                                if (a.length() > 0) {
                                    n++;
                                    name = a.toUpperCase();
                                    if (n > 1) {
                                        break;
                                    }
                                }
                            }
                        }

                        if (readingExons == 1) {
                            int x = line.indexOf(")");
                            if (x > 0) {
                                exonsBuilder.append(line, 0, x);
                                readingExons = 2;
                                System.out.println("Exons: " + exonsBuilder);
                            } else {
                                exonsBuilder.append(line);
                            }
                        }
                        if (readingExons == 0 && line.contains("panel") && line.contains("join(")) {
                            int x = line.indexOf("join(");
                            if (x > 0) {
                                exonsBuilder.append(line.substring(x + 5));
                            }
                            readingExons = 1;
                            x = line.indexOf(")");
                            if (x > 0) {
                                readingExons = 0;
                            }
                        }

                        if (isOriginSection) {
                            sequenceBuilder.append(line);
                        }
                        if (line.contains("origin")) {
                            isOriginSection = true;
                        }
                    }

                    String exons = exonsBuilder.toString();
                    String sequence = sequenceBuilder.toString();

                    if (!sequence.isEmpty() && !exons.isEmpty()) {
                        String[] p = exons.split(",");
                        int[] v = new int[2 * p.length];
                        int n = -1;
                        for (String p1 : p) {
                            String[] r = p1.split("\\.");
                            if (r.length > 1) {
                                v[++n] = StrToInt(r[0]);
                                v[++n] = StrToInt(r[2]);
                            }
                        }
                        String[] newlistprimers = Run(tagfile, sequence, name, primers, v, minpcr, maxpcr, minlen, maxlen, mintm, maxtm, minlc, prlap, ftail, rtail, e5, e3);
                        // Combine the primers with new set of primers
                        primers = combineArrays(primers, newlistprimers);
                    }
                } catch (IOException e) {
                }
            }
        }
    }

    private static String[] Run(String tagfile, String seq, String name, String[] listprimers, int[] exons, int minpcr, int maxpcr, int minlen, int maxlen, int mintm, int maxtm, int minlc, int prlap, String ftail, String rtail, String e5, String e3) {
        String primerlistfile = tagfile + "_primers.txt";
        String pcrcolfile = tagfile + "_panels.txt";
        String panelprimerlistfile = tagfile + "_panelprimers.txt";
        String[] listprimers2 = null;

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
            sr.append("min length=").append(maxlen).append("\n");
            sr.append("Linguistic Complexity(%)=").append(minlc).append("\n");

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

            PrimerDesign pd = new PrimerDesign(dna.DNA(seq), name, minlen, maxlen, mintm, maxtm, minlc, exons, listprimers);
            List<Pair> pcrcol = pd.RunDesign(ftail, rtail, e5, e3, prlap, minpcr, maxpcr);

            PanelsCollector panels = new PanelsCollector(pcrcol);
            panels.CombinePairsPanel(exons.length);
            List<Pair> pcrpanel1 = panels.getPanel1();
            List<Pair> pcrpanel2 = panels.getPanel2();
            ArrayList<String> stringList = new ArrayList<>();

            if (!pcrpanel1.isEmpty()) {
                sr1.append("Panel1:\n");
                for (Pair pcrx : pcrpanel1) {
                    stringList.add(pcrx.fprimer);
                    stringList.add(pcrx.rprimer);
                    srpanel.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\n");
                    srpanel.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\n");
                    sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
                    sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.rlc).append("\n");
                    sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp\n\n");
                }
            }
            if (!pcrpanel2.isEmpty()) {
                sr1.append("Panel2:\n");
                for (Pair pcrx : pcrpanel2) {
                    stringList.add(pcrx.fprimer);
                    stringList.add(pcrx.rprimer);
                    srpanel.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\n");
                    srpanel.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\n");
                    sr1.append(pcrx.fprimername).append("\t").append(pcrx.fprimer).append("\t").append(pcrx.fln).append("\t").append(String.format("%.1f", pcrx.fTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.flc).append("\n");
                    sr1.append(pcrx.rprimername).append("\t").append(pcrx.rprimer).append("\t").append(pcrx.rln).append("\t").append(String.format("%.1f", pcrx.rTm)).append("\t").append(String.format("%.1f", pcrx.fCG)).append("\t").append(pcrx.rlc).append("\n");
                    sr1.append("ExonID:").append((1 + pcrx.exonid)).append(" PCR amplicon=").append(pcrx.pcrsize).append(" bp\n\n");
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
            for (int j = 0; j < exons.length - 1; j += 2) {
                h++;
                sr.append("\n").append("exon:").append(h + 1).append(" ").append(exons[j]).append("-").append(exons[j + 1]).append(" ").append(exons[j + 1] - exons[j] + 1).append("bp\n");
                System.out.println("\nexon:" + (h + 1) + " " + exons[j] + "-" + exons[j + 1] + " " + (exons[j + 1] - exons[j] + 1) + "bp");

                PrimersCollector fPrimersList1 = fPrimersList[h];
                if (fPrimersList1.Amount() > 0) {
                    StringBuilder sr2 = new StringBuilder(100000);
                    sr2.append("Forward:\n");
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
                    sr2.append("\n");
                    System.out.println(sr2);
                    sr.append(sr2);
                }

                PrimersCollector rPrimersList1 = rPrimersList[h];
                if (fPrimersList1.Amount() > 0) {
                    StringBuilder sr2 = new StringBuilder(100000);
                    sr2.append("Reverse:\n");
                    double[] Tm = rPrimersList1.getTms();
                    double[] CG = rPrimersList1.getCGs();
                    int[] ln = rPrimersList1.getPrimerLengths();          // length
                    int[] lc = rPrimersList1.getPrimerLC();               // Linguistic_Complexity
                    //  int[] x1 = rPrimersList1.getPrimerLocations();    // location x1
                    String[] primer = rPrimersList1.getPrimer();          // primer sequence
                    String[] primername = rPrimersList1.getpPrimerName(); //name
                    for (int i = 0; i < rPrimersList1.Amount(); i++) {
                        sr2.append(primername[i]).append("\t").append(primer[i]).append("\t").append(ln[i]).append("\t").append(String.format("%.1f", Tm[i])).append("\t").append(String.format("%.1f", CG[i])).append("\t").append(lc[i]).append("\n");
                    }
                    sr2.append("\n");
                    System.out.println(sr2);
                    sr.append(sr2);
                }

            }

            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds\n\n");
            try (FileWriter fileWriter = new FileWriter(primerlistfile)) {
                System.out.println("Saving the primer list report to a file: " + primerlistfile);
                fileWriter.write(sr.toString());
                fileWriter.write("Time taken: " + duration + " seconds\n\n");
            }

            if (!pcrpanel1.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(panelprimerlistfile)) {
                    System.out.println("Saving panels PCR primers combinations report to a file: " + panelprimerlistfile);
                    fileWriter.write(srpanel.toString());
                }
            }

            if (!pcrcol.isEmpty()) {
                try (FileWriter fileWriter = new FileWriter(pcrcolfile)) {
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
}
