import java.util.ArrayList;
import java.util.List;

public class PanelsCollector {

    private List<Pair> pcrcol;
    private List<Pair> p1;
    private List<Pair> p2;

    /**
     *
     * @param pcrcol
     */
    public PanelsCollector(List<Pair> pcrcol) {
        this.pcrcol = pcrcol;
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

    public void CombinePairsPanel(int exons) {
    p1 = new ArrayList<>();
    p2 = new ArrayList<>();

    // Ожидаем, что pcrcol отсортирован: сначала по exonid, затем по fx5 (возрастание)
    // Если не так — отсортируйте заранее:
    // pcrcol.sort(Comparator.comparingInt((Pair p) -> p.exonid).thenComparingInt(p -> p.fx5));

    for (int j = 0; j < exons; j++) {
        // 1) находим диапазон пар для текущего экзона j
        int n1 = -1, n2 = -1;
        for (int i = 0; i < pcrcol.size(); i++) {
            Pair pcr = pcrcol.get(i);
            if (pcr.exonid == j) {
                if (n1 == -1) n1 = i;
                n2 = i;
            }
            if (pcr.exonid > j) break;
        }
        if (n1 == -1) continue; // для этого экзона пар нет

    
        List<String> panel1Primers = new ArrayList<>();
        List<String> panel2Primers = new ArrayList<>();
        int panel1_lastRx5 = Integer.MIN_VALUE; // px1
        int panel2_lastRx5 = Integer.MIN_VALUE; // px2
        int panel1_guardF5  = Integer.MIN_VALUE; // pf — динамический порог для f5 в панели 2

        for (int i = n1; i <= n2; i++) {
            Pair pcr = pcrcol.get(i);
            final int f5 = pcr.fx5;  // 5' forward
            final int r5 = pcr.rx5;  // 5' reverse

            boolean placed = false;

            // ---- Панель 1 ----
            if (f5 > panel1_lastRx5) {
                if (passesDimerCheck(pcr, panel1Primers)) {
                    // Добавляем праймеры в «белый список» панели 1
                    panel1Primers.add(pcr.fprimer);
                    panel1Primers.add(pcr.rprimer);

                    panel1_lastRx5 = r5;

                    // вычисляем «сторожевой» порог для панели 2 (как у вас: четверть расстояния)
                    panel1_guardF5 = f5 + (r5 - f5) / 4;

                    pcr.SetPanelID(1);
                    p1.add(pcr);
                    placed = true;
                }
            }

            // ---- Панель 2 ----
            if (!placed) {
                // Ваша логика: f5 > px2 && f5 > pf
                if (f5 > panel2_lastRx5 && f5 > panel1_guardF5) {
                    if (passesDimerCheck(pcr, panel2Primers)) {
                        panel2Primers.add(pcr.fprimer);
                        panel2Primers.add(pcr.rprimer);

                        panel2_lastRx5 = r5;

                        pcr.SetPanelID(2);
                        p2.add(pcr);
                        placed = true;
                    }
                }
            }

            // если не поместилось никуда — просто пропускаем
        }
    }
}

private static boolean passesDimerCheck(Pair pcr, List<String> panelPrimers) {
    if (panelPrimers.isEmpty()) return true;
    // Предполагается, что Oligo.DimersCheck(...) возвращает 0 == "нет димеров"
    if (Oligo.DimersCheck(pcr.fprimer, panelPrimers.toArray(String[]::new)) != 0) return false;
    return Oligo.DimersCheck(pcr.rprimer, panelPrimers.toArray(String[]::new)) == 0;
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
