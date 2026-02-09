import java.util.*;

public class ExonsDesigner {

    private final int minpcr;
    private final int maxpcr;

    public ExonsDesigner(int minpcr, int maxpcr) {
        this.minpcr = minpcr;
        this.maxpcr = maxpcr;
    }

    public final class Exon {
        public final int start;
        public final int end;
        public final int idx; // 0-based

        Exon(int start, int end, int idx) {
            this.start = start;
            this.end = end;
            this.idx = idx;
        }
    }

    public final class PcrFragment {
        public final int start;
        public final int end;
        public final List<Integer> exonIndices;

        PcrFragment(int start, int end, List<Integer> exonIndices) {
            this.start = start;
            this.end = end;
            this.exonIndices = Collections.unmodifiableList(new ArrayList<>(exonIndices));
        }

        public int length() { return end - start; }

        @Override public String toString() {
            return "PCR[" + start + "-" + end + " len=" + length() + " exons=" + exonIndices + "]";
        }
    }

    /**
     * Формирует список ПЦР-фрагментов, не изменяя исходные координаты экзонов.
     * exonsArray: [s1,e1,s2,e2,...]  
     */
    public List<PcrFragment> buildPcrFragments(int[] exonsArray) {
        if (exonsArray == null || exonsArray.length % 2 != 0) {
            throw new IllegalArgumentException("exonsArray must be non-null and even-length [start,end,...]");
        }

        // 1) Переносим экзоны в объекты (исходный массив не модифицируем)
        List<Exon> exons = new ArrayList<>(exonsArray.length / 2);
        for (int i = 0, k = 0; i < exonsArray.length; i += 2, k++) {
            int s = exonsArray[i], e = exonsArray[i + 1];
            if (e <= 0) continue;                    // как и раньше пропуск "выключенных"
            if (s > e) { int t = s; s = e; e = t; }  // нормализация на всякий случай
            exons.add(new Exon(s, e, k));
        }

   //     exons.sort(Comparator.comparingInt(x -> x.start));

        List<PcrFragment> fragments = new ArrayList<>();

        // 2) «Скользящие» фрагменты: для каждого i расширяем вправо по правилам
        for (int i = 0; i < exons.size(); i++) {
            int fragStart = exons.get(i).start;
            int fragEnd   = exons.get(i).end;
            List<Integer> included = new ArrayList<>();
            included.add(exons.get(i).idx);

            for (int j = i + 1; j < exons.size(); j++) {
                Exon next = exons.get(j);
                int proposedEnd = Math.max(fragEnd, next.end);
                int newLen = proposedEnd - fragStart;
                int gap = next.start - fragEnd; // может быть ≤0 при перекрытии

                boolean withinMax = newLen <= maxpcr;
                boolean gapTooSmall = gap < minpcr;

                if (withinMax || gapTooSmall) {
                    fragEnd = proposedEnd;
                    included.add(next.idx);
                } else {
                    break;
                }
            }
            fragments.add(new PcrFragment(fragStart, fragEnd, included));
        }

        return fragments;
    }


}
