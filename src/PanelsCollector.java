import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class PanelsCollector {

    private List<Pair> pcrcol;
    private List<Pair> p1;
    private List<Pair> p2;

    private List<Pair> candidates;                      // pcrcol plus the alternative pairs, for closing gaps
    private List<List<Pair>> pools = new ArrayList<>(); // amplicons placed in each tube, all regions and sequences together
    private List<List<Pair>> sortedPools = new ArrayList<>(); // the same, in report order
    private int mainPools;                              // pools tiled as usual: 2 in multiplex, 1 in single-plex
    private boolean extraTubes;                         // multiplex: what conflicts with both pools goes to Panel3, Panel4, ...
    private boolean allowOverlap;                       // single-plex: amplicons are amplified separately
    private List<Pair> forced = new ArrayList<>();      // pairs placed despite a conflict
    private Map<Integer, List<Integer>> noPairs = new HashMap<>(); // per sequence: targets (1-based) that no candidate pair reaches
    private List<RegionNames> names = new ArrayList<>(); // per sequence: names of its regions in the notes

    /**
     *
     * @param pcrcol
     */
    public PanelsCollector(List<Pair> pcrcol) {
        this(pcrcol, new ArrayList<>());
    }

    /**
     * @param pcrcol      candidate pairs, tiled region by region
     * @param alternative further valid pairs (may repeat pcrcol), used only to
     *                    close what the tiling leaves uncovered
     */
    public PanelsCollector(List<Pair> pcrcol, List<Pair> alternative) {
        this.pcrcol = pcrcol;
        this.candidates = new ArrayList<>(pcrcol);
        Set<String> seen = new HashSet<>();
        for (Pair pcr : pcrcol) {
            seen.add(pcr.seq + ":" + pcr.exonid + ":" + pcr.fx5 + ":" + pcr.rx5);
        }
        for (Pair pcr : alternative) {
            if (seen.add(pcr.seq + ":" + pcr.exonid + ":" + pcr.fx5 + ":" + pcr.rx5)) {
                candidates.add(pcr);
            }
        }
    }

    /**
     *
     * @return
     */
    public List<Pair> getPanel1() {
        return p1;
    }

    /**
     *
     * @return
     */
    public List<Pair> getPanel2() {
        return p2;
    }

    /**
     * Number of tubes: the pools, plus in multiplex the extra tubes (Panel3,
     * Panel4, ...) opened for amplicons that conflict with both pools.
     */
    public int getPanelCount() {
        return sortedPools.size();
    }

    /**
     * Amplicons of tube {@code n}, 1-based (Panel1, Panel2, Panel3, ...).
     */
    public List<Pair> getPanel(int n) {
        return sortedPools.get(n - 1);
    }

    /**
     * Pairs placed despite a conflict in their pool; the conflict is in
     * {@link Pair#note}.
     */
    public List<Pair> getForced() {
        return forced;
    }

    /**
     * Names of the regions of each target sequence (index = Pair.seq), used in
     * the notes; without them a region is named by its number.
     */
    public void setRegionNames(List<RegionNames> names) {
        this.names = names;
    }

    private RegionNames namesOf(int seq) {
        return (seq < names.size() && names.get(seq) != null) ? names.get(seq) : new RegionNames(null, null);
    }

    public void CombinePairsPanel(int[] targets, int npools) {
        CombinePairsPanel(List.of(targets), npools);
    }

    /**
     * Assembles the panel so that no target is dropped while a primer pair for
     * it exists.
     *
     * 1. Regions are tiled in order of increasing number of candidate pairs:
     *    the regions with the fewest alternatives seed the pools, and the
     *    regions with many alternatives are fitted around them afterwards.
     *    With several target sequences this order runs across all of them.
     * 2. Primer dimers and amplicon overlaps are checked against everything
     *    already in the pool, across regions and sequences, because a pool is
     *    one tube.
     * 3. Every target still not covered in full is filled from the unused
     *    pairs, alternatives included: first only where they fit without a
     *    conflict, then by moving the conflicting amplicons to another pool.
     * 4. In multiplex, what still conflicts with both pools goes to extra
     *    tubes (Panel3, Panel4, ...): into the first one where it fits, or
     *    into a new one, so that no tube carries a conflict. In single-plex
     *    such a pair is placed in the one pool with the fewest conflicts, and
     *    flagged.
     *
     * In single-plex every amplicon is amplified on its own, so amplicons may
     * overlap; only primer dimers count as conflicts.
     *
     * @param targets per target sequence (index = Pair.seq), its exon
     *                boundaries before merging, [s1,e1,s2,e2,...]
     * @param npools  2 for the two overlapping multiplex pools, 1 for single-plex
     */
    public void CombinePairsPanel(List<int[]> targets, int npools) {
        mainPools = Math.max(1, npools);
        pools = new ArrayList<>();
        for (int k = 0; k < mainPools; k++) {
            pools.add(new ArrayList<>());
        }
        allowOverlap = mainPools == 1;
        extraTubes = mainPools > 1;
        forced = new ArrayList<>();
        noPairs = new HashMap<>();
        for (Pair pcr : candidates) {
            pcr.SetPanelID(0);
            pcr.note = null;
        }

        // 1) candidate pairs of each region, in tiling order
        TreeMap<Long, List<Pair>> byRegion = new TreeMap<>();
        for (Pair pcr : pcrcol) {
            byRegion.computeIfAbsent(((long) pcr.seq << 32) | pcr.exonid, k -> new ArrayList<>()).add(pcr);
        }
        List<List<Pair>> regions = new ArrayList<>(byRegion.values());
        for (List<Pair> region : regions) {
            region.sort(Comparator.comparingInt(p -> p.fx5));
        }

        // 2) problematic regions first (stable sort: ties stay in exon order)
        regions.sort(Comparator.comparingInt(List::size));
        for (List<Pair> region : regions) {
            tileRegion(region);
        }

        // 3) coverage guarantee, again for the targets with the fewest pairs first
        List<int[]> tv = new ArrayList<>(); // {seq, start, end}
        List<Integer> reaching = new ArrayList<>();
        for (int s = 0; s < targets.size(); s++) {
            int[] tb = targets.get(s);
            int n = 0;
            for (int i = 0; i + 1 < tb.length; i += 2) {
                int[] t = CoverageReport.targetInterval(tb[i], tb[i + 1]);
                if (t == null) {
                    continue;
                }
                n++;
                int r = 0;
                for (Pair pcr : candidates) {
                    if (pcr.seq == s && pcr.fx5 < t[1] && pcr.rx5 + 1 > t[0]) {
                        r++;
                    }
                }
                if (r == 0) {
                    noPairs.computeIfAbsent(s, k -> new ArrayList<>()).add(n);
                }
                tv.add(new int[]{s, t[0], t[1]});
                reaching.add(r);
            }
        }
        List<Integer> order = new ArrayList<>();
        for (int n = 0; n < tv.size(); n++) {
            order.add(n);
        }
        order.sort(Comparator.comparingInt(reaching::get));
        List<int[]> byDifficulty = new ArrayList<>();
        for (int n : order) {
            byDifficulty.add(tv.get(n));
        }
        fillGaps(byDifficulty, FIT);
        fillGaps(byDifficulty, RELOCATE);
        fillGaps(byDifficulty, extraTubes ? EXTRA_TUBE : FORCE);

        sortedPools = new ArrayList<>();
        for (int k = 0; k < pools.size(); k++) {
            sortedPools.add(sortedPool(k));
        }
        p1 = sortedPools.get(0);
        p2 = (mainPools > 1) ? sortedPools.get(1) : new ArrayList<>();
    }

    /**
     * Tiles one region: panel 1 takes non-overlapping amplicons, panel 2 the
     * ones that bridge them.
     */
    private void tileRegion(List<Pair> region) {
        int panel1_lastRx5 = Integer.MIN_VALUE;
        int panel2_lastRx5 = Integer.MIN_VALUE;
        int panel1_guardF5 = Integer.MIN_VALUE; // panel 2 starts past the first quarter of the last panel 1 amplicon

        for (Pair pcr : region) {
            final int f5 = pcr.fx5;  // 5' forward
            final int r5 = pcr.rx5;  // 5' reverse

            if (f5 > panel1_lastRx5 && fits(pcr, 0)) {
                place(pcr, 0);
                panel1_lastRx5 = r5;
                panel1_guardF5 = f5 + (r5 - f5) / 4;
            } else if (mainPools > 1 && f5 > panel2_lastRx5 && f5 > panel1_guardF5 && fits(pcr, 1)) {
                place(pcr, 1);
                panel2_lastRx5 = r5;
            }
        }
    }

    private static final int FIT = 0;        // main pools, only where the pair fits as it is
    private static final int RELOCATE = 1;   // main pools, also by moving the amplicons in its way to another pool
    private static final int EXTRA_TUBE = 2; // multiplex: the first extra tube where the pair fits, else a new tube
    private static final int FORCE = 3;      // single-plex: the pool, also despite a conflict

    /**
     * Closes what is left uncovered in the targets ({seq, start, end}) with
     * pairs not placed yet. In mode FIT a pair goes only into a main pool where
     * it neither overlaps an amplicon nor forms a primer dimer; RELOCATE also
     * moves the amplicons in its way to another main pool where they fit;
     * EXTRA_TUBE puts the pair into the first extra tube where it fits, or
     * opens a new tube for it, preferring pairs that need no new tube; FORCE,
     * when nothing else covers the gap, places the pair with the fewest
     * conflicts (overlaps first, then dimers) and flags it in {@link Pair#note}.
     */
    private void fillGaps(List<int[]> targets, int mode) {
        for (int[] t : targets) {
            while (true) {
                List<int[]> covered = coveredIntervals(t[0]);
                if (CoverageReport.overlap(covered, t[1], t[2]) >= t[2] - t[1]) {
                    break;
                }
                Pair best = null;
                int bestPool = -1;
                int bestGain = 0;
                int[] bestScore = null;        // {overlaps, dimers, amplicons moved, new tube}
                List<Pair> bestBlocking = null;
                int[] bestMoves = null;        // new pool of each blocking amplicon
                for (Pair pcr : candidates) {
                    if (pcr.panelID != 0) {
                        continue;
                    }
                    int gain = gain(pcr, t, covered);
                    if (gain <= 0) {
                        continue;
                    }
                    if (mode == EXTRA_TUBE) {
                        int k = mainPools;
                        while (k < pools.size() && !fits(pcr, k)) {
                            k++;
                        }
                        int[] score = {0, 0, 0, (k < pools.size()) ? 0 : 1};
                        if (best == null || better(score, gain, k, bestScore, bestGain, bestPool)) {
                            best = pcr;
                            bestPool = k;
                            bestGain = gain;
                            bestScore = score;
                        }
                        continue;
                    }
                    for (int k = 0; k < mainPools; k++) {
                        int[] score;
                        List<Pair> blocking = null;
                        int[] moves = null;
                        if (fits(pcr, k)) {
                            score = new int[]{0, 0, 0, 0};
                        } else if (mode == RELOCATE) {
                            blocking = blocking(pcr, k);
                            moves = relocation(blocking, k);
                            if (moves == null) {
                                continue;
                            }
                            score = new int[]{0, 0, blocking.size(), 0};
                        } else if (mode == FORCE) {
                            score = conflicts(pcr, k);
                        } else {
                            continue;
                        }
                        if (best == null || better(score, gain, k, bestScore, bestGain, bestPool)) {
                            best = pcr;
                            bestPool = k;
                            bestGain = gain;
                            bestScore = score;
                            bestBlocking = blocking;
                            bestMoves = moves;
                        }
                    }
                }
                if (best == null) {
                    break; // no pair reaches the rest of this target
                }
                if (bestPool == pools.size()) {
                    pools.add(new ArrayList<>()); // a new extra tube
                }
                if (bestMoves != null) {
                    for (int i = 0; i < bestBlocking.size(); i++) {
                        pools.get(bestPool).remove(bestBlocking.get(i));
                        place(bestBlocking.get(i), bestMoves[i]);
                    }
                } else if (bestScore[0] > 0 || bestScore[1] > 0) {
                    best.note = describeConflicts(best, bestPool);
                    forced.add(best);
                }
                place(best, bestPool);
            }
        }
    }

    /**
     * Fewer overlaps, then fewer dimers, then fewer amplicons moved, then no
     * new tube, then more of the gap covered; finally the smaller main pool,
     * or the first extra tube.
     */
    private boolean better(int[] score, int gain, int pool, int[] bestScore, int bestGain, int bestPool) {
        for (int i = 0; i < score.length; i++) {
            if (score[i] != bestScore[i]) {
                return score[i] < bestScore[i];
            }
        }
        if (gain != bestGain) {
            return gain > bestGain;
        }
        if (pool >= mainPools || bestPool >= mainPools) {
            return pool < bestPool;
        }
        return pools.get(pool).size() < pools.get(bestPool).size();
    }

    /**
     * Amplicons of the pool that overlap the pair or form a dimer with it.
     */
    private List<Pair> blocking(Pair pcr, int pool) {
        List<Pair> b = new ArrayList<>();
        for (Pair q : pools.get(pool)) {
            if (overlaps(pcr, q) || dimer(pcr, q)) {
                b.add(q);
            }
        }
        return b;
    }

    /**
     * Another main pool for each blocking amplicon, where it fits without a
     * conflict, or null when one of them has nowhere to go. Forced amplicons
     * stay where they are.
     */
    private int[] relocation(List<Pair> blocking, int from) {
        if (blocking.isEmpty() || mainPools < 2) {
            return null;
        }
        int[] to = new int[blocking.size()];
        for (int i = 0; i < blocking.size(); i++) {
            Pair y = blocking.get(i);
            if (y.note != null) {
                return null;
            }
            to[i] = -1;
            for (int k = 0; k < mainPools && to[i] < 0; k++) {
                if (k == from || !fits(y, k)) {
                    continue;
                }
                boolean ok = true;
                for (int j = 0; j < i && ok; j++) {
                    ok = to[j] != k || !(overlaps(y, blocking.get(j)) || dimer(y, blocking.get(j)));
                }
                if (ok) {
                    to[i] = k;
                }
            }
            if (to[i] < 0) {
                return null;
            }
        }
        return to;
    }

    private void place(Pair pcr, int pool) {
        pcr.SetPanelID(pool + 1);
        pools.get(pool).add(pcr);
    }

    /**
     * True when the pair neither overlaps an amplicon of the pool nor forms a
     * dimer with any of its primers.
     */
    private boolean fits(Pair pcr, int pool) {
        for (Pair q : pools.get(pool)) {
            if (overlaps(pcr, q)) {
                return false;
            }
        }
        for (Pair q : pools.get(pool)) {
            if (dimer(pcr, q)) {
                return false;
            }
        }
        return true;
    }

    /**
     * {overlapping amplicons, amplicons with a primer dimer, 0, 0} in the pool.
     */
    private int[] conflicts(Pair pcr, int pool) {
        int[] c = new int[4];
        for (Pair q : pools.get(pool)) {
            if (overlaps(pcr, q)) {
                c[0]++;
            }
            if (dimer(pcr, q)) {
                c[1]++;
            }
        }
        return c;
    }

    private String describeConflicts(Pair pcr, int pool) {
        List<String> dimers = new ArrayList<>();
        List<String> overlaps = new ArrayList<>();
        for (Pair q : pools.get(pool)) {
            if (overlaps(pcr, q)) {
                overlaps.add(namesOf(q.seq).design(q.exonid) + " " + (q.fx5 + 1) + "-" + (q.rx5 + 1));
            }
            String[][] combos = {
                {pcr.fprimer, pcr.fprimername, q.fprimer, q.fprimername},
                {pcr.fprimer, pcr.fprimername, q.rprimer, q.rprimername},
                {pcr.rprimer, pcr.rprimername, q.fprimer, q.fprimername},
                {pcr.rprimer, pcr.rprimername, q.rprimer, q.rprimername}};
            for (String[] x : combos) {
                if (dimer(x[0], x[2])) {
                    dimers.add(x[1] + " x " + x[3]);
                }
            }
        }
        StringBuilder s = new StringBuilder();
        if (!dimers.isEmpty()) {
            s.append("primer-dimer ").append(String.join(", ", dimers));
        }
        if (!overlaps.isEmpty()) {
            if (s.length() > 0) {
                s.append("; ");
            }
            s.append("overlaps ").append(String.join(", ", overlaps));
        }
        return s.toString();
    }

    private boolean overlaps(Pair a, Pair b) {
        return !allowOverlap && a.seq == b.seq && a.fx5 <= b.rx5 && b.fx5 <= a.rx5;
    }

    private boolean dimer(Pair a, Pair b) {
        return dimer(a.fprimer, b.fprimer) || dimer(a.fprimer, b.rprimer)
                || dimer(a.rprimer, b.fprimer) || dimer(a.rprimer, b.rprimer);
    }

    // The same primers are checked against each other many times over
    private final Map<String, Integer> primerIds = new HashMap<>();
    private final LongIntHashMap dimerCache = new LongIntHashMap(1 << 16);

    /**
     * Oligo.quickDimer can answer differently with its arguments swapped
     * (primers of equal length), so both orders are checked.
     */
    private boolean dimer(String a, String b) {
        long ia = primerIds.computeIfAbsent(a, k -> primerIds.size());
        long ib = primerIds.computeIfAbsent(b, k -> primerIds.size());
        // LongIntHashMap hashes by folding the two halves of the key, so the pair
        // of ids is scrambled first (SplitMix64 finaliser, a bijection)
        long key = (ia <= ib) ? (ia << 32) | ib : (ib << 32) | ia;
        key = (key ^ (key >>> 30)) * 0xbf58476d1ce4e5b9L;
        key = (key ^ (key >>> 27)) * 0x94d049bb133111ebL;
        key ^= key >>> 31;
        boolean cacheable = key > Long.MIN_VALUE + 1; // the two smallest keys are the map's empty/deleted marks
        int d = cacheable ? dimerCache.get(key) : LongIntHashMap.NO_VALUE;
        if (d == LongIntHashMap.NO_VALUE) {
            d = (Oligo.quickDimer(a, b) != 0 || Oligo.quickDimer(b, a) != 0) ? 1 : 0;
            if (cacheable) {
                dimerCache.put(key, d);
            }
        }
        return d == 1;
    }

    /**
     * Uncovered bases of target {@code t} ({seq, start, end}) that the pair
     * would cover.
     */
    private static int gain(Pair pcr, int[] t, List<int[]> covered) {
        if (pcr.seq != t[0]) {
            return 0;
        }
        int lo = Math.max(t[1], pcr.fx5);
        int hi = Math.min(t[2], pcr.rx5 + 1);
        if (hi <= lo) {
            return 0;
        }
        return (hi - lo) - CoverageReport.overlap(covered, lo, hi);
    }

    private List<int[]> coveredIntervals(int seq) {
        List<int[]> iv = new ArrayList<>();
        for (List<Pair> pool : pools) {
            for (Pair q : pool) {
                if (q.seq == seq) {
                    iv.add(CoverageReport.interval(q.fx5, 1 + q.rx5));
                }
            }
        }
        return CoverageReport.merge(iv);
    }

    private List<Pair> sortedPool(int pool) {
        List<Pair> list = new ArrayList<>(pools.get(pool));
        list.sort(Comparator.comparingInt((Pair p) -> p.seq).thenComparingInt(p -> p.exonid).thenComparingInt(p -> p.fx5));
        return list;
    }

    /**
     * Report section for one target sequence on what the coverage guarantee
     * had to do; empty when nothing needed an extra tube or was forced, and
     * every target has a candidate pair.
     */
    public String renderNotes(int seq) {
        List<Pair> seqForced = new ArrayList<>();
        for (Pair pcr : forced) {
            if (pcr.seq == seq) {
                seqForced.add(pcr);
            }
        }
        List<String> extra = new ArrayList<>(); // "Panel3 (5)": the extra tubes this sequence uses
        for (int k = mainPools; k < pools.size(); k++) {
            int n = 0;
            for (Pair pcr : pools.get(k)) {
                if (pcr.seq == seq) {
                    n++;
                }
            }
            if (n > 0) {
                extra.add("Panel" + (k + 1) + " (" + n + ")");
            }
        }
        List<Integer> seqNoPairs = noPairs.getOrDefault(seq, new ArrayList<>());
        if (extra.isEmpty() && seqForced.isEmpty() && seqNoPairs.isEmpty()) {
            return "";
        }
        StringBuilder s = new StringBuilder();
        s.append("\n\nPanel assembly\n");
        if (!extra.isEmpty()) {
            s.append("Amplicons that conflict with both Panel1 and Panel2 are amplified in extra tubes, ")
                    .append("each free of conflicts, so that their regions are not dropped: ")
                    .append(String.join(", ", extra)).append("\n");
            s.append("The panel as a whole needs ").append(pools.size()).append(" tubes.\n");
        }
        if (!seqForced.isEmpty()) {
            s.append("Amplicons placed despite a conflict, so that their region is not dropped: ")
                    .append(seqForced.size()).append("\n");
            s.append("Check them before ordering.\n");
            for (Pair pcr : seqForced) {
                s.append(namesOf(seq).design(pcr.exonid)).append("\tPanel").append(pcr.panelID)
                        .append("\t").append(pcr.fprimername).append("\t").append(pcr.rprimername)
                        .append("\t").append(pcr.note).append("\n");
            }
        }
        if (!seqNoPairs.isEmpty()) {
            s.append("Regions without any candidate primer pair at the current settings: ");
            for (int i = 0; i < seqNoPairs.size(); i++) {
                if (i > 0) {
                    s.append(", ");
                }
                s.append(namesOf(seq).target(seqNoPairs.get(i) - 1));
            }
            s.append("\n");
        }
        return s.toString();
    }

    public void CombinePairsPanel2(int exons) {
        p1 = new ArrayList<>();
        p2 = new ArrayList<>();
        String[] pr1 = new String[0]; // list primers pane 1
        String[] pr2 = new String[0]; // list primers pane 2
        int px1 = 0;       // x5 -reverse primer panel 1 
        int px2 = 0;       // x5 -reverse primer panel 2 
        int pf = 0;// x5 -forward primer panel 1 
        int pr = 0;// x5 -reverse primer panel 1         

        for (int j = 0; j < exons; j++) {
            Pair pcr;
            int n1 = -1;
            int n2 = -1;
            for (int i = 0; i < pcrcol.size(); i++) {
                pcr = pcrcol.get(i);
                if (pcr.exonid == j) {
                    if (n1 == -1) {
                        n1 = i;
                    }
                    n2 = i;
                }
                if (pcr.exonid > j) {
                    break;
                }
            }

            if (n1 > -1 && n2 > -1) {

                for (int i = n1; i <= n2; i++) {
                    pcr = pcrcol.get(i);
                    int f5 = pcr.fx5;
                    int q1 = 1;

                    if (f5 > px1) //adding to panel 1
                    {
                        if (pr1.length > 0) {
                            q1 = Oligo.DimersCheck(pcr.fprimer, pr1);
                            if (q1 == 0) {
                                q1 = Oligo.DimersCheck(pcr.rprimer, pr1);
                            }
                            if (q1 == 0) {
                                int k = pr1.length;
                                String[] s = new String[k + 2];
                                System.arraycopy(pr1, 0, s, 0, k);
                                pr1 = s;
                                pr1[k] = pcr.fprimer;
                                pr1[k + 1] = pcr.rprimer;
                                px1 = pcr.rx5;

                                pr = pcr.rx5;
                                pf = f5 + (pr - f5) / 4;

                                pcr.SetPanelID(1);
                                p1.add(pcr);
                            }
                        } else {
                            pr1 = new String[2];
                            pr1[0] = pcr.fprimer;
                            pr1[1] = pcr.rprimer;
                            px1 = pcr.rx5;

                            pr = pcr.rx5;
                            pf = f5 + (pr - f5) / 4;

                            pcr.SetPanelID(1);
                            p1.add(pcr);
                            q1 = 0;

                        }
                    }

                    if (q1 > 0) { // panel 2
                        if (f5 > px2 && f5 > pf) {//&& f5 > pf) {
                            if (pr2.length > 0) {
                                q1 = Oligo.DimersCheck(pcr.fprimer, pr2);
                                if (q1 == 0) {
                                    q1 = Oligo.DimersCheck(pcr.rprimer, pr2);
                                }
                                if (q1 == 0) {
                                    int k = pr2.length;
                                    String[] s = new String[k + 2];
                                    System.arraycopy(pr2, 0, s, 0, k);
                                    pr2 = s;
                                    pr2[k] = pcr.fprimer;
                                    pr2[k + 1] = pcr.rprimer;
                                    px2 = pcr.rx5;
                                    pcr.SetPanelID(2);
                                    p2.add(pcr);
                                }
                            } else {
                                pr2 = new String[2];
                                pr2[0] = pcr.fprimer;
                                pr2[1] = pcr.rprimer;
                                px2 = pcr.rx5;
                                pcr.SetPanelID(2);
                                p2.add(pcr);
                            }
                        }
                    }
                }
            }
        }
    }
}
