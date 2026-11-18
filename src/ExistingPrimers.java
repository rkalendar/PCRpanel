import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * The primers a panel already has, read from the {@code target_primers} list.
 *
 * Two primers on consecutive lines that face each other on a target sequence
 * are an existing pair: its amplicon counts as covering what it reads, and the
 * panel is designed only for what such pairs leave unread. This is the layout
 * of {@code *_panelprimers.txt}, a forward primer followed by its reverse one.
 * A primer may carry an adapter tail; its 3' part is what has to match.
 *
 * A line holds a name and a sequence, separated by a tab or a space, or a
 * sequence alone; a third column may give the tube (1, 2, ...) the primer is
 * amplified in. Primers that are in no pair on any target stay what the list
 * was before: oligos that no new primer may form a dimer with.
 */
public class ExistingPrimers {

    private static final int MIN_MATCH = 15;      // bases at the 3' end that have to match the target
    private static final int MAX_PRODUCT = 10000; // bp: two primers further apart than this are not a pair

    private static final class Entry {

        String name;
        String seq;
        int tube;       // 0 = not given
        boolean paired; // in a pair on some target sequence
    }

    private final List<Entry> entries = new ArrayList<>();
    private final primer calculator = new primer();
    private final List<String> tails = new ArrayList<>();

    /**
     * The adapter tails of the run. A primer that begins with one binds with
     * the rest of it; without this, the last bases of a tail that happen to
     * match the target would be taken for part of the binding site.
     */
    public void setTails(String... tails) {
        this.tails.clear();
        for (String t : tails) {
            if (t != null && !t.isEmpty()) {
                this.tails.add(t.toLowerCase());
            }
        }
    }

    private String withoutTail(String primer) {
        for (String t : tails) {
            if (primer.startsWith(t) && primer.length() - t.length() >= MIN_MATCH) {
                return primer.substring(t.length());
            }
        }
        return primer;
    }

    /** Reads the list; an empty or missing file name gives an empty list. */
    public static ExistingPrimers read(String file) {
        ExistingPrimers list = new ExistingPrimers();
        if (file == null || file.length() <= 1) {
            return list;
        }
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                // the sequence is the first column made of bases only; what precedes it is the name
                String[] f = line.split("[\\t ]+");
                int at = (f.length == 1) ? 0 : -1;
                for (int k = 1; k < f.length && at < 0; k++) {
                    if (f[k].matches("(?i)[acgtunrymkswbdhvi]{9,}")) {
                        at = k;
                    }
                }
                if (at < 0) {
                    continue;
                }
                Entry e = new Entry();
                e.name = (at > 0) ? String.join("", java.util.Arrays.copyOfRange(f, 0, at)) : "primer" + (list.entries.size() + 1);
                e.seq = dna.DNA(f[at].toLowerCase());
                if (at + 1 < f.length && f[at + 1].matches("\\d+")) {
                    e.tube = Integer.parseInt(f[at + 1]);
                }
                if (e.seq.length() > 8) {
                    list.entries.add(e);
                }
            }
        } catch (IOException e) {
            System.err.println("Warning: could not read primer file '" + file + "': " + e.getMessage());
        }
        return list;
    }

    public boolean isEmpty() {
        return entries.isEmpty();
    }

    public int size() {
        return entries.size();
    }

    /**
     * The pairs of the list that lie on {@code sequence}: consecutive lines, in
     * either order, whose 3' parts match the two strands facing each other.
     */
    public List<Pair> pairsOn(String sequence) {
        List<Pair> pairs = new ArrayList<>();
        if (entries.size() < 2) {
            return pairs;
        }
        String antisense = dna.ComplementDNA(sequence);
        int n = sequence.length();
        List<List<int[]>> forward = new ArrayList<>(); // binding sites on the sense strand, {start, end}
        List<List<int[]>> reverse = new ArrayList<>();
        for (Entry e : entries) {
            String binding = withoutTail(e.seq);
            forward.add(sites(binding, sequence, false, n));
            reverse.add(sites(binding, antisense, true, n));
        }
        for (int i = 0; i + 1 < entries.size(); i++) {
            int[][] fr = facing(forward.get(i), reverse.get(i + 1));
            int f = i;
            int r = i + 1;
            if (fr == null) {
                fr = facing(forward.get(i + 1), reverse.get(i));
                f = i + 1;
                r = i;
            }
            if (fr == null) {
                continue;
            }
            Entry a = entries.get(f);
            Entry b = entries.get(r);
            Pair p = new Pair(0, 0, -1, fr[1][1] - fr[0][0]);
            String fb = a.seq.substring(a.seq.length() - (fr[0][1] - fr[0][0]));
            p.AddFprimer(a.seq, a.name, fr[0][0], fb.length(), calculator.getLC(fb), calculator.getTm(fb), calculator.getGC());
            String rb = b.seq.substring(b.seq.length() - (fr[1][1] - fr[1][0]));
            p.AddRprimer(b.seq, b.name, fr[1][1] - 1, rb.length(), calculator.getLC(rb), calculator.getTm(rb), calculator.getGC());
            p.existing = true;
            p.SetPanelID((a.tube == b.tube) ? a.tube : 0);
            a.paired = true;
            b.paired = true;
            pairs.add(p);
            i++; // both lines are used
        }
        return pairs;
    }

    /**
     * Where the primer binds a strand: every site of the longest 3' part of it,
     * at least MIN_MATCH long, that the strand holds. A tail at its 5' end
     * need not match.
     */
    private static List<int[]> sites(String primer, String strand, boolean antisense, int n) {
        List<int[]> out = new ArrayList<>();
        for (int len = primer.length(); len >= MIN_MATCH && out.isEmpty(); len--) {
            String part = primer.substring(primer.length() - len);
            for (int at = strand.indexOf(part); at >= 0; at = strand.indexOf(part, at + 1)) {
                out.add(antisense ? new int[]{n - at - len, n - at} : new int[]{at, at + len});
            }
        }
        return out;
    }

    /** The forward and the reverse site giving the shortest product, or null if none face each other. */
    private static int[][] facing(List<int[]> forward, List<int[]> reverse) {
        int[][] best = null;
        for (int[] f : forward) {
            for (int[] r : reverse) {
                int size = r[1] - f[0];
                if (f[1] <= r[0] && size <= MAX_PRODUCT && (best == null || size < best[1][1] - best[0][0])) {
                    best = new int[][]{f, r};
                }
            }
        }
        return best;
    }

    /** The primers that are in no pair on any target sequence so far. */
    public String[] unpaired() {
        List<String> out = new ArrayList<>();
        for (Entry e : entries) {
            if (!e.paired) {
                out.add(e.seq);
            }
        }
        return out.toArray(new String[0]);
    }

    /** Their names, for the console. */
    public List<String> unpairedNames() {
        List<String> out = new ArrayList<>();
        for (Entry e : entries) {
            if (!e.paired) {
                out.add(e.name);
            }
        }
        return out;
    }
}
