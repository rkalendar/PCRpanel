import java.util.*;

public class PCRcollector {

    List<Map.Entry<Integer, Integer>> entryList = new ArrayList<>();
    PrimersCollector[] fPrimersList;
    PrimersCollector[] rPrimersList;
    int currentSize;
    int exon1;
    int exon2;

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
     * @return all candidate pairs, region by region
     */
    public List<Pair> CombinePrimers(int minpcr, int maxpcr) {
        List<Pair> c = new ArrayList<>();
        for (int i = 0; i < currentSize; i++) {
            c.addAll(pairPrimers(fPrimersList[i], rPrimersList[i], minpcr, maxpcr));
        }
        return c;
    }

    /**
     * Every valid combination of forward and reverse primers, region by region:
     * the alternatives the panel can fall back on to close a coverage gap.
     */
    public List<Pair> CombineAllPrimers(int minpcr, int maxpcr) {
        List<Pair> c = new ArrayList<>();
        for (int i = 0; i < currentSize; i++) {
            c.addAll(pairAllPrimers(fPrimersList[i], rPrimersList[i], minpcr, maxpcr));
        }
        return c;
    }

    /**
     * All pairs of one region that {@link #pairPrimers} could form, without
     * limiting a forward or reverse primer to a single pair.
     */
    public static List<Pair> pairAllPrimers(PrimersCollector fp, PrimersCollector rp, int minpcr, int maxpcr) {
        List<Pair> c = new ArrayList<>();
        if (fp == null || rp == null || fp.Amount() == 0 || rp.Amount() == 0) {
            return c;
        }
        int exonid = (fp.getTaskId() - 1);
        int exon1 = fp.getTaskX1();
        int exon2 = fp.getTaskX2();
        if (exon2 - exon1 < minpcr) {
            exon1 = fp.getTaskX2();
            exon2 = fp.getTaskX1();
        }
        int[] fx1 = fp.getPrimerLocations();
        int[] rx1 = rp.getPrimerLocations();
        String[] fprimer = fp.getPrimer();
        String[] rprimer = rp.getPrimer();

        for (int f = 0; f < fprimer.length; f++) {
            for (int r = 0; r < rprimer.length; r++) {
                int tpcr = rx1[r] - fx1[f];
                if (rx1[r] > 0 && rx1[r] > exon1 && fx1[f] < exon2 && tpcr >= minpcr && tpcr <= maxpcr
                        && Oligo.quickDimer(fprimer[f], rprimer[r]) == 0) {
                    Pair p = new Pair(f, r, exonid, 1 + tpcr);
                    p.AddFprimer(fprimer[f], fp.getpPrimerName()[f], fx1[f], fp.getPrimerLengths()[f], fp.getPrimerLC()[f], fp.getTms()[f], fp.getCGs()[f]);
                    p.AddRprimer(rprimer[r], rp.getpPrimerName()[r], rx1[r], rp.getPrimerLengths()[r], rp.getPrimerLC()[r], rp.getTms()[r], rp.getCGs()[r]);
                    c.add(p);
                }
            }
        }
        return c;
    }

    /**
     * Candidate pairs of one region: each forward primer takes the nearest unused
     * reverse primer that gives an amplicon of minpcr..maxpcr overlapping the
     * region (spanning it, for a region shorter than minpcr).
     *
     * Primers of other regions are deliberately not checked here. Checking every
     * region against all candidates of the regions before it starved the later
     * regions of pairs, although most of those candidates never reach the panel;
     * PanelsCollector checks dimers against the primers actually placed in a pool.
     */
    public static List<Pair> pairPrimers(PrimersCollector fp, PrimersCollector rp, int minpcr, int maxpcr) {
        List<Pair> c = new ArrayList<>();
        if (fp == null || rp == null || fp.Amount() == 0 || rp.Amount() == 0) {
            return c;
        }
        int exonid = (fp.getTaskId() - 1);
        int exon1 = fp.getTaskX1();          // exon x1
        int exon2 = fp.getTaskX2();          // exon x2
        if (exon2 - exon1 < minpcr) {
            exon1 = fp.getTaskX2();
            exon2 = fp.getTaskX1();
        }

        double[] fTm = fp.getTms();
        double[] fCG = fp.getCGs();
        int[] fln = fp.getPrimerLengths();          // length
        int[] flc = fp.getPrimerLC();               // Linguistic_Complexity
        int[] fx1 = fp.getPrimerLocations();        // location x1
        String[] fprimer = fp.getPrimer();          // primer sequence
        String[] fprimername = fp.getpPrimerName(); // name

        double[] rTm = rp.getTms();
        double[] rCG = rp.getCGs();
        int[] rln = rp.getPrimerLengths();          // length
        int[] rlc = rp.getPrimerLC();               // Linguistic_Complexity
        int[] rx1 = rp.getPrimerLocations().clone(); // location x1; a used reverse primer is set to -1
        String[] rprimer = rp.getPrimer();          // primer sequence
        String[] rprimername = rp.getpPrimerName(); // name

        for (int f = 0; f < fprimer.length; f++) {
            int x5f = fx1[f]; // 5-end forward
            for (int r = rprimer.length - 1; r > -1; r--) { // reverse primers in increasing position
                int x5r = rx1[r]; // 5-end reverse
                if (x5r - x5f > maxpcr) {
                    break;
                }
                if (x5r > 0 && x5r > exon1 && x5f < exon2) {
                    int tpcr = x5r - x5f;
                    if (tpcr >= minpcr && tpcr <= maxpcr && Oligo.quickDimer(fprimer[f], rprimer[r]) == 0) {
                        Pair p = new Pair(f, r, exonid, 1 + tpcr);//Pair(int forward, int reverse, int exonid, int pcrsize) {
                        p.AddFprimer(fprimer[f], fprimername[f], fx1[f], fln[f], flc[f], fTm[f], fCG[f]);
                        p.AddRprimer(rprimer[r], rprimername[r], rx1[r], rln[r], rlc[r], rTm[r], rCG[r]);
                        c.add(p);
                        rx1[r] = -1;
                        break;
                    }
                }
            }
        }
        return c;
    }
}
