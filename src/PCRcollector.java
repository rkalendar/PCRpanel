
import java.util.*;

public class PCRcollector {

    List<Map.Entry<Integer, Integer>> entryList = new ArrayList<>();
    PrimersCollector[] fPrimersList;
    PrimersCollector[] rPrimersList;
    int currentSize;
    int exon1;
    int exon2;
    int minamplicon = 50;

    public PCRcollector() {

        currentSize = 0;
        fPrimersList = new PrimersCollector[1];
        rPrimersList = new PrimersCollector[1];
    }

    public void add(PrimersCollector f, PrimersCollector r) {
        currentSize++;
        PrimersCollector[] newFPrimersList = new PrimersCollector[currentSize];
        PrimersCollector[] newRPrimersList = new PrimersCollector[currentSize];
        System.arraycopy(fPrimersList, 0, newFPrimersList, 0, fPrimersList.length);
        System.arraycopy(rPrimersList, 0, newRPrimersList, 0, rPrimersList.length);
        fPrimersList = newFPrimersList;
        rPrimersList = newRPrimersList;
        fPrimersList[currentSize - 1] = f;
        rPrimersList[currentSize - 1] = r;
        entryList.add(new AbstractMap.SimpleEntry<>(currentSize - 1, Math.min(f.Amount(), r.Amount())));
    }

    public int getCurrentSize() {
        entryList.sort(Map.Entry.comparingByValue());
        return currentSize;
    }

    /**
     *
     * @param minpcr
     * @param maxpcr
     * @return
     */
    public List<Pair> CombinePrimers(int minpcr, int maxpcr) {
        Set<String> totalset = new HashSet<>();
        List<Pair> c = new ArrayList<>();
        String[] primerslist = new String[0];
        List<Integer> zs = new ArrayList<>();

        for (int i = 0; i < currentSize; i++) {
            int exonid = (fPrimersList[i].getTaskId() - 1);
            int exon1 = fPrimersList[i].getTaskX1();          // exon x1
            int exon2 = fPrimersList[i].getTaskX2();          // exon x2
            if (exon2 - exon1 < minpcr) {
                exon1 = fPrimersList[i].getTaskX2();
                exon2 = fPrimersList[i].getTaskX1();
            }

            double[] fTm = fPrimersList[i].getTms();
            double[] fCG = fPrimersList[i].getCGs();
            int[] fln = fPrimersList[i].getPrimerLengths();          // length
            int[] flc = fPrimersList[i].getPrimerLC();               // Linguistic_Complexity
            int[] fx1 = fPrimersList[i].getPrimerLocations();        // location x1
            String[] fprimer = fPrimersList[i].getPrimer();          // primer sequence
            String[] fprimername = fPrimersList[i].getpPrimerName(); // name

            double[] rTm = rPrimersList[i].getTms();
            double[] rCG = rPrimersList[i].getCGs();
            int[] rln = rPrimersList[i].getPrimerLengths();          // length
            int[] rlc = rPrimersList[i].getPrimerLC();               // Linguistic_Complexity
            int[] rx1 = rPrimersList[i].getPrimerLocations();        // location x1
            String[] rprimer = rPrimersList[i].getPrimer();          // primer sequence
            String[] rprimername = rPrimersList[i].getpPrimerName(); // name

            for (int f = 0; f < fprimer.length; f++) {
                int x5f = fx1[f]; // 5-end forward
                for (int r = rprimer.length - 1; r > -1; r--) { // for (int r = 0; r <rprimer.length ; r++) { 
                    int x5r = rx1[r]; // 5-end reverse
                    if (1 + x5r - x5f >= maxpcr) {
                        zs.add(1 + x5r - x5f);
                        break;
                    }
                    if (x5r - x5f > minamplicon) {
                        if (x5r > exon1 && x5f < exon2) {
                            zs.add(1 + x5r - x5f);
                        }
                    }
                }
            }

            int mnpcr = minpcr;
            int mxpcr = maxpcr;
            if (!zs.isEmpty()) {
                Collections.sort(zs);
                if (zs.getFirst() < minpcr && zs.getLast() > maxpcr) {
                    //mnpcr = zs.getFirst();
                    mxpcr = zs.getLast();
                }
            }

            for (int f = 0; f < fprimer.length; f++) {
                int x5f = fx1[f]; // 5-end forward
                for (int r = rprimer.length - 1; r > -1; r--) { // for (int r = 0; r <rprimer.length ; r++) { 
                    int x5r = rx1[r]; // 5-end reverse
                    if (x5r - x5f > mxpcr) {
                        break;
                    }

                    if (rx1[r] > 0) {
                        if (x5r > exon1 && x5f < exon2) {
                            int tpcr = x5r - x5f;
                            if (tpcr >= mnpcr & tpcr <= maxpcr) {
                                if (Oligo.quickDimer(fprimer[f], rprimer[r]) == 0) {
                                    int y = 0;
                                    if (primerslist.length > 0 && Oligo.DimersCheck(fprimer[f], primerslist) > 0) {
                                        y = 1;
                                    }
                                    if (y == 0) {
                                        Pair p = new Pair(f, r, exonid, 1 + tpcr);//Pair(int forward, int reverse, int exonid, int pcrsize) {
                                        p.AddFprimer(fprimer[f], fprimername[f], fx1[f], fln[f], flc[f], fTm[f], fCG[f]);
                                        p.AddRprimer(rprimer[r], rprimername[r], rx1[r], rln[r], rlc[r], rTm[r], rCG[r]);
                                        c.add(p);
                                        totalset.add(fprimer[f]);
                                        totalset.add(rprimer[r]);
                                        rx1[r] = -1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (!totalset.isEmpty()) {
                primerslist = totalset.toArray(String[]::new);
            }
        }
        return c;
    }
}
