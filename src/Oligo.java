import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public final class Oligo {

    public static double e260(String str) {
        int l = str.length();
        if (l < 4) {
            return (0);
        }
        int n1 = str.charAt(0);
        int n2 = str.charAt(1);
        double t = oligoparam.vODdna[n1 - 97][n2 - 97];
        for (int i = 1; i < l - 1; i++) {
            n1 = str.charAt(i);
            n2 = str.charAt(i + 1);
            t = t + oligoparam.vODdna[n1 - 97][n2 - 97] + oligoparam.nODdna[n1 - 97];
        }
        return (t);
    }

    public static int Gdimer(String r1, String r2, String m1[], String sb[], String m2[]) {
        int i = 0;
        int f = 0;
        int j = 0;
        int w = 0;
        int z = 0;
        int n = 0;
        int x = 0;
        int x1 = 0;
        int x2 = 0;
        int n1 = 0;
        int n2 = 0;
        String z2 = r1;

        int l1 = r1.length();
        int l2 = r2.length();

        if (l2 > l1) {
            r1 = r2;
            r2 = z2;
            l1 = r1.length();
            l2 = r2.length();
        }

        String[] s1 = {"ggg", "ccc"};
        char[] s2 = {'g', 'c'};

        int[][] p1 = new int[l1 - 2][2];
        int[][] p2 = new int[l2 - 2][2];

        String c = tools.Strings(l1 + l2, ' ');
        z2 = dna.ReverseSeq(r2);

        for (int r = 0; r < 2; r++) {
            for (i = 0;; i++) {
                i = r1.indexOf(s1[r], i);
                if (i == -1) {
                    break;
                }
                n1++;
                p1[n1][0] = i;
                p1[n1][1] = 3;
                j = i + 3;
                for (i = j; i < l1 && r1.charAt(i) == s2[r]; i++) {
                    p1[n1][1]++;
                }
            }

            for (i = -1;; i++) {
                i = z2.indexOf(s1[r], i);
                if (i == -1) {
                    break;
                }
                n2++;
                p2[n2][0] = i;
                p2[n2][1] = 3;
                j = i + 3;
                for (i = j; i < l2 && z2.charAt(i) == s2[r]; i++) {
                    p2[n2][1]++;
                }
            }

            if (n1 > 1 & n2 > 1) {
                for (i = 1; i < n1; i++) {
                    if (p1[i + 1][0] - p1[i][0] - p1[i][1] < 2 + Math.min(p1[i + 1][1], p1[i][1])) {
                        for (j = 1; j < n2; j++) {
                            if (p2[j + 1][0] - p2[j][0] - p2[j][1] < 2 + Math.min(p2[j + 1][1], p2[j][1])) {

                                for (n = 0; n < Math.min(Math.min(p1[i][1] - 3, p2[j][1] - 3), Math.min(p1[i + 1][1] - 3, p2[j + 1][1] - 3)); n++) {
                                    if (p2[j][1] > p1[i][1]) {
                                        x1 = p1[i + 1][0] - p1[i][0] - p1[i][1] - n;
                                        x2 = p2[j + 1][0] - p2[j][0] - p2[j][1];
                                        w = p1[i][0] - p2[j][0] - n;
                                    } else {
                                        x1 = p1[i + 1][0] - p1[i][0] - p1[i][1];
                                        x2 = p2[j + 1][0] - p2[j][0] - p2[j][1] - n;
                                        w = p1[i][0] - p2[j][0] + n;
                                    }
                                    if (x1 > x2) {
                                        z = x2 + p2[j + 1][1] - x1;
                                    } else {
                                        z = x1 + p1[i + 1][1] - x2;
                                    }
                                    if (z > 2) {
                                        f++;
                                        m1[f] = "<-" + z2 + "-5";
                                        m2[f] = "5-" + r1 + "->";
                                        if (w > 0) {
                                            m1[f] = c.substring(0, w) + "<-" + z2 + "-5";
                                            m2[f] = "5-" + r1 + "->";
                                        }
                                        if (w < 0) {
                                            m1[f] = "<-" + z2 + "-5";
                                            m2[f] = c.substring(0, -w) + "5-" + r1 + "->";
                                        }
                                        int l = Math.min(m1[f].length(), m2[f].length());
                                        StringBuilder value = new StringBuilder(c.substring(0, l));
                                        for (x = 2; x < l - 2; x++) {
                                            char a1 = m1[f].charAt(x);
                                            char a2 = m2[f].charAt(x);
                                            int h = tables.aDNA[a1][a2];
                                            if (h > 70) {
                                                value.setCharAt(x, '|');
                                            } else {
                                                if (h > tables.smmin) {
                                                    value.setCharAt(x, ':');
                                                }
                                                if (a1 == a2 && (tables.dx[a1] > 1) && (m1[f].charAt(x - 1) == m2[f].charAt(x - 1) || m1[f].charAt(x + 1) == m2[f].charAt(x + 1))) {
                                                    value.setCharAt(x, '|');
                                                }
                                            }
                                        }
                                        sb[f] = new String(value);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return f;
    }

    public static int Gdimer(String r1, String r2) {
        int i = -1;
        int f = 0;
        int j = 0;
        int z = 0;
        int n = 0;
        int x1 = 0;
        int x2 = 0;

        int l1 = r1.length();
        int l2 = r2.length();

        String z2 = dna.ReverseSeq(r2);
        String[] s1 = {"ggg", "ccc"};
        char[] s2 = {'g', 'c'};

        int[][] p1 = new int[l1 - 2][2];
        for (int r = 0; r < 2; r++) {
            int n1 = 0;
            for (i = 0;; i++) {
                i = r1.indexOf(s1[r], i);
                if (i == -1) {
                    break;
                }
                n1++;
                p1[n1][0] = i;
                p1[n1][1] = 3;
                j = i + 3;
                for (i = j; i < l1 && r1.charAt(i) == s2[r]; i++) {
                    p1[n1][1]++;
                }
            }
            int n2 = 0;
            int[][] p2 = new int[l2 - 2][2];
            for (i = 0;; i++) {
                i = z2.indexOf(s1[r], i);
                if (i == -1) {
                    break;
                }
                n2++;
                p2[n2][0] = i;
                p2[n2][1] = 3;
                j = i + 3;
                for (i = j; i < l2 && z2.charAt(i) == s2[r]; i++) {
                    p2[n2][1]++;
                }
            }

            if (n1 > 1 & n2 > 1) {
                for (i = 1; i < n1; i++) {
                    if (p1[i + 1][0] - p1[i][0] - p1[i][1] < 2 + Math.min(p1[i + 1][1], p1[i][1])) {
                        for (j = 1; j < n2; j++) {
                            if (p2[j + 1][0] - p2[j][0] - p2[j][1] < 2 + Math.min(p2[j + 1][1], p2[j][1])) {
                                for (n = 0; n < Math.min(Math.min(p1[i][1] - 3, p2[j][1] - 3), Math.min(p1[i + 1][1] - 3, p2[j + 1][1] - 3)); n++) {
                                    if (p2[j][1] > p1[i][1]) {
                                        x1 = p1[i + 1][0] - p1[i][0] - p1[i][1] - n;
                                        x2 = p2[j + 1][0] - p2[j][0] - p2[j][1];
                                    } else {
                                        x1 = p1[i + 1][0] - p1[i][0] - p1[i][1];
                                        x2 = p2[j + 1][0] - p2[j][0] - p2[j][1] - n;
                                    }
                                    if (x1 > x2) {
                                        z = x2 + p2[j + 1][1] - x1;
                                    } else {
                                        z = x1 + p1[i + 1][1] - x2;
                                    }
                                    if (z > 2) {
                                        f++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if ((n1 == 1 & n2 == 1) && (p1[1][1] > 6 && p2[1][1] > 6)) {
                f++;
            }
        }
        return f;
    }

    public static int Gdimer(String r1, String m1[], String sb[], String m2[]) {
        int f = 0;
        int j = 0;
        int w = 0;
        int x = 0;
        int l = r1.length();
        String[] s1 = {"ggg", "ccc"};
        char[] s2 = {'g', 'c'};
        int[][] p1 = new int[l - 2][2];

        String c = tools.Strings(l + l, ' ');
        String z2 = dna.ReverseSeq(r1);

        for (int r = 0; r < 2; r++) {
            int n = 0;
            int i = 0;
            for (;; i++) {
                i = r1.indexOf(s1[r], i);
                if (i == -1) {
                    break;
                }
                n++;
                p1[n][0] = i;
                p1[n][1] = 3;
                j = i + 3;
                for (i = j; i < l && r1.charAt(i) == s2[r]; i++) {
                    p1[n][1]++;
                }
            }

            if (n > 1) {
                for (i = 1; i < n; i++) {
                    if (p1[i + 1][0] - p1[i][0] - p1[i][1] < 2 + Math.min(p1[i + 1][1], p1[i][1])) {
                        w = p1[i][0] - l + p1[i + 1][0] + p1[i + 1][1];
                        f++;
                        m1[f] = "<-" + z2 + "-5";
                        m2[f] = "5-" + r1 + "->";
                        if (w > 0) {
                            m1[f] = c.substring(0, w) + "<-" + z2 + "-5";
                            m2[f] = "5-" + r1 + "->";
                        }
                        if (w < 0) {
                            m1[f] = "<-" + z2 + "-5";
                            m2[f] = c.substring(0, -w) + "5-" + r1 + "->";
                        }
                        int l1 = Math.min(m1[f].length(), m2[f].length());
                        StringBuilder value = new StringBuilder(c.substring(0, l1));
                        for (x = 2; x < l1 - 2; x++) {
                            char a1 = m1[f].charAt(x);
                            char a2 = m2[f].charAt(x);
                            if (tables.aDNA[a1][a2] > 70) {
                                value.setCharAt(x, '|');
                            } else {
                                if (tables.aDNA[a1][a2] > tables.smmin) {
                                    value.setCharAt(x, ':');
                                }
                                if (a1 == a2 && (tables.dx[a1] > 1) && (m1[f].charAt(x - 1) == m2[f].charAt(x - 1) || m1[f].charAt(x + 1) == m2[f].charAt(x + 1))) {
                                    value.setCharAt(x, '|');
                                }
                            }
                        }
                        sb[f] = new String(value);
                    }
                }
            }
            if (n == 1 && p1[1][1] > 6) {
                w = p1[1][0] - (l - (p1[1][0] + p1[1][1]));
                f++;
                m1[f] = "<-" + z2 + "-5";
                m2[f] = "5-" + r1 + "->";
                if (w > 0) {
                    m1[f] = c.substring(0, w) + "<-" + z2 + "-5";
                    m2[f] = "5-" + r1 + "->";
                }
                if (w < 0) {
                    m1[f] = "<-" + z2 + "-5";
                    m2[f] = c.substring(0, -w) + "5-" + r1 + "->";
                }
                int l1 = Math.min(m1[f].length(), m2[f].length());
                StringBuilder value = new StringBuilder(c.substring(0, l1));
                for (x = 2; x < l1 - 2; x++) {
                    char a1 = m1[f].charAt(x);
                    char a2 = m2[f].charAt(x);
                    int h = tables.aDNA[a1][a2];
                    if (h > 70) {
                        value.setCharAt(x, '|');
                    } else {
                        if (h > tables.smmin) {
                            value.setCharAt(x, ':');
                        }
                        if (a1 == a2 && (tables.dx[a1] > 1) && (m1[f].charAt(x - 1) == m2[f].charAt(x - 1) || m1[f].charAt(x + 1) == m2[f].charAt(x + 1))) {
                            value.setCharAt(x, '|');
                        }
                    }
                }
                sb[f] = new String(value);
            }
        }
        return f;
    }

    public static int Gdimer(String r1) {
        int f = 0;
        int j = 0;
        int l = r1.length();
        String[] s1 = {"ggg", "ccc"};
        char[] s2 = {'g', 'c'};
        int[][] p1 = new int[l][2];
        for (int r = 0; r < 2; r++) {
            int n = 0;
            int i = 0;
            for (;; i++) {
                i = r1.indexOf(s1[r], i);
                if (i == -1) {
                    break;
                }
                n++;
                p1[n][0] = i;
                p1[n][1] = 3;
                j = i + 3;
                for (i = j; i < l && r1.charAt(i) == s2[r]; i++) {
                    p1[n][1]++;
                }
            }
            if (n > 1) {
                for (i = 1; i < n; i++) {
                    if (p1[i + 1][0] - p1[i][0] - p1[i][1] < 2 + Math.min(p1[i + 1][1], p1[i][1])) {
                        f++;
                    }
                }
            }
            if (n == 1 && p1[1][1] > 6) {
                f++;
            }
        }
        return f;
    }

    public static int DimerLook(String r1, String r2, int n3, int n5, int m0[], String m1[], String sb[], String m2[]) {
        int f = 0;
        int n = 0;
        int w = 0;
        int x = 0;
        int y = 0;
        int l = 0;
        int h = 0;
        int g = 0;
        int j = 0;
        int e = 0;
        int b = 0;

        int y1 = 0;  // real y    
        int y0 = 0;
        int y2 = 0;
        int fn = 0;
        int x0 = 0;
        int x2 = 0;
        int xk = 0;
        int x1 = 0;  // real x    

        int b1 = 0;
        int b2 = 0;

        if (n5 < 1) {
            n5 = 99999;
        }

        boolean e3;
        boolean e5;
        //   r1 = "tctacagtccgacgatcgccacttcc";
        //   r2 = r1;
        //   r1 = "caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct";
        //   r2 = "gatcggaagagcggttcagcaggaatgccgag";

        //  r1 = "caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgat";
        // r2 = "gatcggaagagcggttcagcaggaatgccgagaaaaa";
        String z1 = "";
        String z2 = "";
        int l1 = r1.length(); // longest
        int l2 = r2.length(); // shortest

        r1 = r1.toLowerCase();
        r2 = r2.toLowerCase();

        int lmin = l2;
        if (lmin > l1) {
            lmin = l1;
            l1 = l2;
            l2 = lmin;
            z1 = r2;
            z2 = dna.ReverseSeq(r1);
        } else {
            z1 = r1;
            z2 = dna.ReverseSeq(r2);
        }
        if (lmin < 6) {
            return (0);
        }

        int[] lx = new int[l1 + l2 + l2];
        int[] u = new int[l1 - 4];
        int[][][][][][] px = new int[5][5][5][5][5][2];
        int[] ax = new int[5];

        ax[0] = tables.dx[z1.charAt(0)];
        ax[1] = tables.dx[z1.charAt(1)];
        ax[2] = tables.dx[z1.charAt(2)];
        ax[3] = tables.dx[z1.charAt(3)];
        for (j = 4; j < l1; j++) {
            ax[4] = tables.dx[z1.charAt(j)];
            px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]++;
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }
        ax[0] = tables.dx[z1.charAt(0)];
        ax[1] = tables.dx[z1.charAt(1)];
        ax[2] = tables.dx[z1.charAt(2)];
        ax[3] = tables.dx[z1.charAt(3)];
        for (j = 4; j < l1; j++) {
            ax[4] = tables.dx[z1.charAt(j)];
            y = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
            if (y > 0) {
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] = fn;
                u[fn] = j - 4;
                fn = fn + y;
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0] = -1;
            } else {
                u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] - y] = j - 4;
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]--;
            }
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }

        ax[0] = tables.dx[tables.cdna[z2.charAt(0)]];
        ax[1] = tables.dx[tables.cdna[z2.charAt(1)]];
        ax[2] = tables.dx[tables.cdna[z2.charAt(2)]];
        ax[3] = tables.dx[tables.cdna[z2.charAt(3)]];
        for (y = 0; y < l2 - 4; y++) {
            ax[4] = tables.dx[tables.cdna[z2.charAt(y + 4)]];
            f = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
            if (f < 0) {
                for (h = 1; h <= -f; h++) {
                    x = u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] + h - 1];
                    y1 = y + 4;
                    x1 = x + 4;
                    w = x - y;              // location 3' and x
                    xk = l2 + w - 1;        // location 5'-l2 at l1
                    b = l1 - l2 - x - y;

                    if (lx[xk] == 0) {
                        e3 = false;
                        e5 = false;

                        if (w < 1) {
                            x0 = 0;
                            y0 = -w;
                        } else {
                            x0 = w;
                            y0 = 0;
                        }
                        if (b < 0) {
                            x2 = l1;
                            y2 = l2 + b;
                        } else {
                            x2 = l1 - b;
                            y2 = l2;
                        }

//'                            <-gagccgtaaggacgacttggcgagaaggctag-5
//'                              ||||||||||||||||||||||||||||||||
//'5-caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct->
                        if ((x0 > 0 && y0 == 0) || (y2 < l2 && x2 == l1)) { // 3'end control
                            l = 0;
                            g = 1;
                            if (w > 0 && y == 0) {
                                for (j = 0; j < l2 - w; j++) {
                                    b1 = z1.charAt(x0 + j);
                                    b2 = z2.charAt(y0 + j);
                                    if (tables.aDNA[b1][b2] > tables.smmax) {
                                        if (tables.dx[b1] > 1) {
                                            g++;
                                        }
                                        l++;
                                    } else {
                                        break;
                                    }
                                }
                                e3 = (l > n3 || g > n3); // g=5
                            }

//'                       <-aaaaagagccgtaaggacgacttggcgagaaggctag-5
//'                              |||||||||||||||||||||||||||||||
//'5-caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgat->                            
                            if (y2 < l2 && x2 == l1 && !e3) {
                                l = 0;
                                g = 1;
                                for (j = 1; j < y2 + b - 1; j++) {
                                    b1 = z1.charAt(x2 - j);
                                    b2 = z2.charAt(y2 - j);
                                    if (tables.aDNA[b1][b2] > tables.smmax) {
                                        if (tables.dx[b1] > 1) {
                                            g++;
                                        }
                                        l++;
                                    } else {
                                        break;
                                    }
                                }
                                e3 = (l > n3 || g >= n3);
                            }
                        }

                        if (!e3) {
                            l = 5;
                            e = 0;
                            g = 0;
                            if (ax[0] > 1) {
                                g++;
                            }
                            if (ax[1] > 1) {
                                g++;
                            }
                            if (ax[2] > 1) {
                                g++;
                            }
                            if (ax[3] > 1) {
                                g++;
                            }
                            if (ax[4] > 1) {
                                g++;
                            }

                            for (j = 1; j < Math.min(l1 - x1, l2 - y1); j++) {
                                b1 = z1.charAt(x1 + j);
                                b2 = z2.charAt(y1 + j);
                                if (tables.aDNA[b1][b2] == 100) {
                                    if (tables.dx[b1] > 1) {
                                        g++;
                                    }
                                    l++;
                                } else {
                                    if (e == 1) {
                                        break;
                                    }
                                    e = 1;
                                    if (tables.aDNA[b1][b2] > tables.smmin) {
                                        l++;
                                    }
                                }
                            }

                            for (j = 0; j < Math.min(x, y); j++) {
                                b1 = z1.charAt(x - j - 1);
                                b2 = z2.charAt(y - j - 1);
                                if (tables.aDNA[b1][b2] == 100) {
                                    l++;
                                    if (tables.dx[b1] > 1) {
                                        g++;
                                    }
                                } else {
                                    if (e == 1) {
                                        break;
                                    }
                                    e = 1;
                                    if (tables.aDNA[b1][b2] > tables.smmin) {
                                        l++;
                                    }
                                }
                            }
                            e5 = (l > n5 || g > n3);
                        }

                        if (e3 || e5) {
                            lx[xk] = 1;
                            n++;
                            m1[n] = "<-" + z2 + "-5";
                            sb[n] = "  ";
                            m2[n] = "5-" + z1 + "->";
                            m0[n] = 2;
                            if (e3) {
                                m0[n] = 1;
                            }
                            if (w > 0) {
                                m1[n] = tools.Strings(w, ' ') + "<-" + z2 + "-5";
                                sb[n] = tools.Strings(w + 2, ' ');
                                m2[n] = "5-" + z1 + "->";
                            }
                            if (w == 0) {
                                m1[n] = "<-" + z2 + "-5";
                                sb[n] = tools.Strings(2, ' ');
                                m2[n] = "5-" + z1 + "->";
                            }
                            if (w < 0) {
                                m1[n] = "<-" + z2 + "-5";
                                sb[n] = tools.Strings(2 - w, ' ');
                                m2[n] = tools.Strings(-w, ' ') + "5-" + z1 + "->";
                            }
                            int z = Math.min(l2 - y0, l1 - x0);
                            char[] value = new char[z];
                            for (j = 0; j < z; j++) {
                                value[j] = ' ';
                                int h1 = tables.aDNA[z2.charAt(y0 + j)][z1.charAt(x0 + j)];
                                if (h1 > tables.smmin) {
                                    value[j] = ':';
                                }
                                if (h1 > 70) {
                                    value[j] = '|';
                                }
                            }
                            sb[n] = sb[n] + new String(value);
                        }
                    }
                }
            }
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }

        String[] g1 = new String[l1 + l2];
        String[] g2 = new String[l1 + l2];
        String[] g3 = new String[l1 + l2];
        w = Gdimer(r1, g1, g2, g3);
        if (w > 0) {
            for (j = 1; j < w + 1; j++) {
                n++;
                m1[n] = g1[j];
                m2[n] = g3[j];
                sb[n] = g2[j];
                m0[n] = 3;
            }
        }
        return n;
    }

    public static int DimerLookList1(String[] r, String t, int n3, int n5) {
        if (t == null || t.isEmpty() || r == null) {
            return 0;
        }
        for (String r1 : r) {
            int k = DimerLook(r1, t, n3, n5);
            if (k > 0) {
                return k;
            }
        }
        return 0;
    }

    public static int DimerLook(String r1, String r2, int n3, int n5) {
        if (r1 == null || r2 == null || r1.isEmpty() || r2.isEmpty()) {
            return 0;
        }
        int f = 0;
        int n = 0;
        int w = 0;
        int x = 0;
        int y = 0;
        int l = 0;
        int h = 0;
        int j = 0;
        int g = 0;
        int e = 0;
        int b = 0;
        int y1 = 0;  // real y    
        int y0 = 0;
        int y2 = 0;
        int fn = 0;
        int x0 = 0;
        int x2 = 0;
        int xk = 0;
        int x1 = 0;  // real x    

        int b1 = 0;
        int b2 = 0;

        String z1 = "";
        String z2 = "";
        int l1 = r1.length(); // longest
        int l2 = r2.length(); // shortest
        if (n5 == -1) {
            n5 = 99999;
        }

        int lmin = l2;
        if (lmin > l1) {
            lmin = l1;
            l1 = l2;
            l2 = lmin;
            z1 = r2;
            z2 = dna.ReverseSeq(r1);
        } else {
            z1 = r1;
            z2 = dna.ReverseSeq(r2);
        }
        if (lmin < 6) {
            return (0);
        }

        int[] lx = new int[l1 + l2 + l2];
        int[] u = new int[l1 - 4];
        int[][][][][][] px = new int[5][5][5][5][5][2];
        int[] ax = new int[5];

        ax[0] = tables.dx[z1.charAt(0)];
        ax[1] = tables.dx[z1.charAt(1)];
        ax[2] = tables.dx[z1.charAt(2)];
        ax[3] = tables.dx[z1.charAt(3)];
        for (j = 4; j < l1; j++) {
            ax[4] = tables.dx[z1.charAt(j)];
            px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]++;
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }
        ax[0] = tables.dx[z1.charAt(0)];
        ax[1] = tables.dx[z1.charAt(1)];
        ax[2] = tables.dx[z1.charAt(2)];
        ax[3] = tables.dx[z1.charAt(3)];
        for (j = 4; j < l1; j++) {
            ax[4] = tables.dx[z1.charAt(j)];
            y = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
            if (y > 0) {
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] = fn;
                u[fn] = j - 4;
                fn = fn + y;
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0] = -1;
            } else {
                u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] - y] = j - 4;
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]--;
            }
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }

        ax[0] = tables.dx[tables.cdna[z2.charAt(0)]];
        ax[1] = tables.dx[tables.cdna[z2.charAt(1)]];
        ax[2] = tables.dx[tables.cdna[z2.charAt(2)]];
        ax[3] = tables.dx[tables.cdna[z2.charAt(3)]];
        for (y = 0; y < l2 - 4; y++) {
            ax[4] = tables.dx[tables.cdna[z2.charAt(y + 4)]];
            f = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
            if (f < 0) {
                for (h = 1; h <= -f; h++) {
                    x = u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] + h - 1];
                    w = x - y;              // location 3' and x
                    xk = l2 + w - 1;        // location 5'-l2 at l1
                    b = l1 - l2 - x - y;

                    if (lx[xk] == 0) {
                        if (w < 1) {
                            x0 = 0;
                            y0 = -w;
                        } else {
                            x0 = w;
                            y0 = 0;
                        }
                        if (b < 0) {
                            x2 = l1;
                            y2 = l2 + b;
                        } else {
                            x2 = l1 - b;
                            y2 = l2;
                        }

//'                            <-gagccgtaaggacgacttggcgagaaggctag-5
//'                              ||||||||||||||||||||||||||||||||
//'5-caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct->
                        if ((x0 > 0 && y0 == 0) || (y2 < l2 && x2 == l1)) { // 3'end control
                            l = 0;
                            g = 1;
                            if (w > 0 || (b == 0 && w == 0)) {
                                for (j = 0; j < l2 - w; j++) {
                                    b1 = z1.charAt(x0 + j);
                                    b2 = z2.charAt(y0 + j);
                                    if (tables.aDNA[b1][b2] > tables.smmax) {
                                        if (tables.dx[b1] > 1) {
                                            g++;
                                        }
                                        l++;
                                    } else {
                                        break;
                                    }
                                }
                                if (l > n3 || g >= n3) {
                                    lx[xk] = 1;
                                    n++;
                                }
                            }

                            if (y2 < l2 && x2 == l1 && lx[xk] == 0) {
                                l = 0;
                                g = 1;
                                for (j = 1; j < y2 + b - 1; j++) {
                                    b1 = z1.charAt(x2 - j);
                                    b2 = z2.charAt(y2 - j);
                                    if (tables.aDNA[b1][b2] > tables.smmax) {
                                        if (tables.dx[b1] > 1) {
                                            g++;
                                        }
                                        l++;
                                    } else {
                                        break;
                                    }
                                }

                                if (l > n3 || g >= n3) {
                                    lx[xk] = 1;
                                    n++;
                                }
                            }

                        }

                        if (lx[xk] == 0) {
                            l = 5;
                            e = 0;
                            g = 0;
                            if (ax[0] > 1) {
                                g++;
                            }
                            if (ax[1] > 1) {
                                g++;
                            }
                            if (ax[2] > 1) {
                                g++;
                            }
                            if (ax[3] > 1) {
                                g++;
                            }
                            if (ax[4] > 1) {
                                g++;
                            }

                            for (j = 1; j < Math.min(l1 - x1, l2 - y1); j++) {
                                b1 = z1.charAt(x1 + j);
                                b2 = z2.charAt(y1 + j);
                                if (tables.aDNA[b1][b2] == 100) {
                                    if (tables.dx[b1] > 1) {
                                        g++;
                                    }
                                    l++;
                                } else {
                                    if (e == 1) {
                                        break;
                                    }
                                    e = 1;
                                    if (tables.aDNA[b1][b2] > tables.smmin) {
                                        l++;
                                    }
                                }
                            }

                            for (j = 0; j < Math.min(x, y); j++) {
                                b1 = z1.charAt(x - j - 1);
                                b2 = z2.charAt(y - j - 1);
                                if (tables.aDNA[b1][b2] == 100) {
                                    l++;
                                    if (tables.dx[b1] > 1) {
                                        g++;
                                    }
                                } else {
                                    if (e == 1) {
                                        break;
                                    }
                                    e = 1;
                                    if (tables.aDNA[b1][b2] > tables.smmin) {
                                        l++;
                                    }
                                }
                            }
                            if (l > n5) {
                                lx[xk] = 1;
                                n++;
                            }

                        }
                    }
                }
            }
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }
        n = n + Gdimer(r1);
        return n;
    }

    public static int quickDimer(String r1) {
        return qDimer(r1, r1);
    }

    public static int quickDimer(String r1, String r2) {
        return qDimer(r1, r2);
    }

    public static int qDimer(String p1, String p2) {
        final int dimer5 = 8;
        final int dimer3 = 7;
        final int n = 5;

        String r1 = p1;
        String r2 = p2;
        if (p1.length() > p2.length()) {
            r1 = p2;
            r2 = p1;
        }
     /*  
        r1 = "gaaccgaagcgtacagtcgcc";
        r2 = "actccaaaggcgacatctcc";
        /*  
          <-acaaacagtttgtattcgtacgt-5 
            |||||| |||||| 
5-tgcatgcttatgtttgacaaaca->                      
         */
 
        Map<String, Set<Integer>> x1 = new HashMap<>();
        Map<String, Set<Integer>> r1Map = new HashMap<>();
        for (int y = 0; y <= r1.length() - n; y++) {
            String sub = r1.substring(y, y + n);
            r1Map.computeIfAbsent(sub, k -> new HashSet<>()).add(y);
        }

         String z2 = dna.ComplementDNA(r2);
        for (int y = 0; y <= r2.length() - n; y++) {
            String sub = z2.substring(y, y + n);
            if (r1Map.containsKey(sub)) {
                Set<Integer> x2 = r1Map.get(sub);
                int[] x = x2.stream().mapToInt(Integer::intValue).toArray();

                for (int i = 0; i < x.length; i++) {
                    int x0 = r1.length() - x[i] - n;
                    int y0 = z2.length() - y - n;
                    int l = n;
                    int e = 0;

                    if (x0 == 1 && y0 > 1) {
                        for (int k = 0; k < Math.min(y, x[i]); k++) {
                            String s1 = r1.substring(x[i] - k, x[i] - k + 1);
                            String s2 = z2.substring(y - k, y - k + 1);
                            if (s1.equals(s2)) {
                                l++;
                            } else {
                                if (++e > 1) {
                                    break;
                                }
                            }
                        }
                        /*  
                 x[i] x0      
5-gaaccgaagcgtacagtcgcc-> r1
            |    ||||||         
         <-cctctacagcggaaacctca-5        
                 y    y0       
                         */
                    } else {
                        for (int k = 0; k < Math.min(x0, y0); k++) {
                            String s1 = r1.substring(x[i] + k + n, x[i] + k + n + 1);
                            String s2 = z2.substring(y + k + n, y + k + n + 1);
                            if (s1.equals(s2)) {
                                l++;
                            } else {
                                if (++e > 1) {
                                    break;
                                }
                            }
                        }
                    }

                    if (l > dimer3) {
                        if (l > dimer5) {   // 5-end dimer
                            return 1;
                        }
                        if (y == 0 && x[i] > 0) { //3-end dimer
                            return 1;
                        }
                        if (x0 == 1 && y0 > 1) { //3-end dimer
                            return 1;
                        }
                    }
                }
            }
        }
        return 0;
    }

    public static int quickDimerLook(String r1, String r2, int n3) {
  /*
         r1 = "caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct";
         r2 = "gatcggaagagcggttcagcaggaatgccgag";
*/
        int n = 0;
        if (n3 < 6) {
            n3 = 6;
        }
        if (n3 > 10) {
            n3 = 10;
        }

        int l1 = r1.length(); // longest
        int l2 = r2.length(); // shortest
        String z1 = "";
        String z2 = "";

        int lmin = l2;
        if (lmin > l1) {
            lmin = l1;
            l1 = l2;
            l2 = lmin;
            z1 = r2;
            z2 = dna.ComplementDNA(r1);
        } else {
            z1 = r1;
            z2 = dna.ComplementDNA(r2);
        }
        if (lmin < 6) {
            return (0);
        }

        HashMap<String, int[]> map = new HashMap<>();

        for (int y = 0; y < r1.length() - n3 + 1; y++) {
            String s = r1.substring(y, y + n3);
            if (map.containsKey(s)) {
                int[] p = map.get(s);
                int[] v = new int[p.length + 1];
                System.arraycopy(p, 0, v, 0, p.length);
                p = v;
                p[p.length - 1] = y;
                map.replace(s, p);
            } else {
                int[] p = new int[1];
                p[0] = y;
                map.put(s, p);
            }
        }

        String s = z1.substring(0, n3);
        if (map.containsKey(s)) {
            int[] p = map.get(s);
            for (int i = 0; i < p.length; i++) {
                int x = p[i];
            }
        }

        return (int) n;
    }

    public static double DimerLookTm(String s, int n3, int n5) {
        return DimerLookTm(s, s, n3, n5, 1);
    }

    public static double DimerLookTm(String r1, String r2, int n3, int n5, int k) {
        if (r1 == null || r2 == null || r1.isEmpty() || r2.isEmpty()) {
            return 0;
        }
        int f = 0;
        double n = 0;
        int w = 0;
        int x = 0;
        int y = 0;
        int l = 0;
        int h = 0;
        int j = 0;
        int g = 0;
        int e = 0;
        int b = 0;
        int y1 = 0;  // real y    
        int y0 = 0;
        int y2 = 0;
        int fn = 0;
        int x0 = 0;
        int x2 = 0;
        int xk = 0;
        int x1 = 0;  // real x    

        int b1 = 0;
        int b2 = 0;
        String z1 = "";
        String z2 = "";
        int l1 = r1.length(); // longest
        int l2 = r2.length(); // shortest
        if (n5 == -1) {
            n5 = 99999;
        }

        int lmin = l2;
        if (lmin > l1) {
            lmin = l1;
            l1 = l2;
            l2 = lmin;
            z1 = r2;
            z2 = dna.ReverseSeq(r1);
        } else {
            z1 = r1;
            z2 = dna.ReverseSeq(r2);
        }
        if (lmin < 6) {
            return (0);
        }
        String c = tools.Strings(l1 + l2, ' ');

        int[] lx = new int[l1 + l2 + l2];
        int[] u = new int[l1 - 4];
        int[][][][][][] px = new int[5][5][5][5][5][2];
        int[] ax = new int[5];

        ax[0] = tables.dx[z1.charAt(0)];
        ax[1] = tables.dx[z1.charAt(1)];
        ax[2] = tables.dx[z1.charAt(2)];
        ax[3] = tables.dx[z1.charAt(3)];
        for (j = 4; j < l1; j++) {
            ax[4] = tables.dx[z1.charAt(j)];
            px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]++;
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }
        ax[0] = tables.dx[z1.charAt(0)];
        ax[1] = tables.dx[z1.charAt(1)];
        ax[2] = tables.dx[z1.charAt(2)];
        ax[3] = tables.dx[z1.charAt(3)];
        for (j = 4; j < l1; j++) {
            ax[4] = tables.dx[z1.charAt(j)];
            y = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
            if (y > 0) {
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] = fn;
                u[fn] = j - 4;
                fn = fn + y;
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0] = -1;
            } else {
                u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] - y] = j - 4;
                px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]--;
            }
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }

        ax[0] = tables.dx[tables.cdna[z2.charAt(0)]];
        ax[1] = tables.dx[tables.cdna[z2.charAt(1)]];
        ax[2] = tables.dx[tables.cdna[z2.charAt(2)]];
        ax[3] = tables.dx[tables.cdna[z2.charAt(3)]];
        for (y = 0; y < l2 - 4; y++) {
            ax[4] = tables.dx[tables.cdna[z2.charAt(y + 4)]];
            f = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
            if (f < 0) {
                for (h = 1; h <= -f; h++) {
                    x = u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] + h - 1];
                    w = x - y;              // location 3' and x
                    xk = l2 + w - 1;        // location 5'-l2 at l1
                    b = l1 - l2 - x - y;

                    if (lx[xk] == 0) {
                        if (w < 1) {
                            x0 = 0;
                            y0 = -w;
                        } else {
                            x0 = w;
                            y0 = 0;
                        }
                        if (b < 0) {
                            x2 = l1;
                            y2 = l2 + b;
                        } else {
                            x2 = l1 - b;
                            y2 = l2;
                        }

//'                            <-gagccgtaaggacgacttggcgagaaggctag-5
//'                              ||||||||||||||||||||||||||||||||
//'5-caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct->
                        if ((x0 > 0 && y0 == 0) || (y2 < l2 && x2 == l1)) { // 3'end control
                            l = 0;
                            g = 0;
                            if (w > 0 || (b == 0 && w == 0)) {
                                for (j = 0; j < l2 - w; j++) {
                                    b1 = z1.charAt(x0 + j);
                                    b2 = z2.charAt(y0 + j);
                                    if (tables.aDNA[b1][b2] > tables.smmax) {
                                        if (tables.dx[b1] > 1) {
                                            g++;
                                        }
                                        l++;
                                    } else {
                                        break;
                                    }
                                }
                                if (l > n3 || g >= n3) {
                                    lx[xk] = 1;
                                }
                            }

//'                       <-aaaaagagccgtaaggacgacttggcgagaaggctag-5
//'                              |||||||||||||||||||||||||||||||
//'5-caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgat->                            
                            if (y2 < l2 && x2 == l1 && lx[xk] == 0) {
                                for (j = 1; j < y2 + b - 1; j++) {
                                    b1 = z1.charAt(x2 - j);
                                    b2 = z2.charAt(y2 - j);
                                    if (tables.aDNA[b1][b2] > tables.smmax) {
                                        if (tables.dx[b1] > 1) {
                                            g++;
                                        }
                                        l++;
                                    } else {
                                        break;
                                    }
                                }
                            }
                            if (l > n3 || g >= n3) {
                                lx[xk] = 1;

                                String s1 = z2;
                                String s2 = z1;
                                if (w > 0) {
                                    s1 = c.substring(0, w) + z2;
                                }
                                if (w < 0) {
                                    s2 = c.substring(0, -w) + z1;
                                }
                                double tm = Melting.Tmelting(s2, s1, 0.25d, 50d, 0, 0, k);
                                if (tm > n) {
                                    n = tm;
                                }
                            }
                        }// 3-end 

                        if (lx[xk] == 0) {
                            l = 5;
                            e = 0;
                            g = 0;
                            if (ax[0] > 1) {
                                g++;
                            }
                            if (ax[1] > 1) {
                                g++;
                            }
                            if (ax[2] > 1) {
                                g++;
                            }
                            if (ax[3] > 1) {
                                g++;
                            }
                            if (ax[4] > 1) {
                                g++;
                            }

                            for (j = 1; j < Math.min(l1 - x1, l2 - y1); j++) {
                                b1 = z1.charAt(x1 + j);
                                b2 = z2.charAt(y1 + j);
                                if (tables.aDNA[b1][b2] == 100) {
                                    if (tables.dx[b1] > 1) {
                                        g++;
                                    }
                                    l++;
                                } else {
                                    if (e == 1) {
                                        break;
                                    }
                                    e = 1;
                                    if (tables.aDNA[b1][b2] > tables.smmin) {
                                        l++;
                                    }
                                }
                            }

                            for (j = 0; j < Math.min(x, y); j++) {
                                b1 = z1.charAt(x - j - 1);
                                b2 = z2.charAt(y - j - 1);
                                if (tables.aDNA[b1][b2] == 100) {
                                    l++;
                                    if (tables.dx[b1] > 1) {
                                        g++;
                                    }
                                } else {
                                    if (e == 1) {
                                        break;
                                    }
                                    e = 1;
                                    if (tables.aDNA[b1][b2] > tables.smmin) {
                                        l++;
                                    }
                                }
                            }
                            if (l > n5) {
                                lx[xk] = 1;

                                String s1 = z2;
                                String s2 = z1;
                                if (w > 0) {
                                    s1 = c.substring(0, w) + z2;
                                }
                                if (w < 0) {
                                    s2 = c.substring(0, -w) + z1;
                                }
                                double tm = Melting.Tmelting(s2, s1, 0.25d, 50d, 0, 0, k);
                                if (tm > n) {
                                    n = tm;
                                }
                            }

                        }
                    }
                }
            }
            ax[0] = ax[1];
            ax[1] = ax[2];
            ax[2] = ax[3];
            ax[3] = ax[4];
        }

        return n;
    }

    public static int PrimerQuality(String source, int x) {
        return PrimerQuality(source, x, false);
    }

    public static int PrimerQuality(String source) {
        return PrimerQuality(source, 0, false);
    }

    public static int PrimerQuality(String source, int dm, boolean probe) {
        int l = source.length();
        double e3 = 0;
        double e5 = 0;
        double r = 0;

        primer prm = new primer(source, 55, 1, 0, 0.2);// (String primer, double Ksalt_mM, double mgcl_mM, double dmso0, double primerconc) {
        double n = prm.getDegeneracy();
        if (n > 0) {
            n = 7 * Math.log(n);
        }
        n = n + (5 * prm.CountTaq("ggg"));
        n = n + (5 * prm.CountTaq("cccc"));
        n = n + (5 * prm.CountTaq("aaaa"));
        n = n + (5 * prm.CountTaq("tttt"));
        if (dm == 0) {
            n = n + DimerLookTm(source, tables.min3dimer, tables.min5dimer);
        }

        if (probe) {
            /* 1. S1 left-hand sequence from the target, in a temperature range below 40;
   2. S2 Right side of the probe, temperature below 30;
   3. Tm for probe as a whole to be about 55.             */
            int x = 12;
            if (l > 12) {
                x = (l * 2) / 3;
            }
            String s = source.substring(0, x);
            e5 = prm.getTm(s);
            s = source.substring(l - x);
            e3 = prm.getTm(s);
            r = e5 + e3 - 10;
            if (e3 + 5 < e5) {
                r = (e5 - e3) * 5;
            }

        } else {
            String s = source.substring(l - 12);
            e3 = prm.getTm(s);
            s = source.substring(0, 12);
            e5 = prm.getTm(s);
            r = (e3 < 51 ? e3 : 50 + (2 * (50 - e3)));
            r = r + r;
            r = r + (e5 < 44 ? e5 : 87 - e5);
            r = r - 40;
        }

        r = r - n;
        if (r < 0) {
            r = 0;
        }
        return ((int) r);
    }

    

    public static int Ta(int l, double t1, double t2) {
        if (l < 50) {
            l = 50;
        }
        double t = (t2 < t1) ? t2 : t1;
        int r = (int) (t + Math.log(l));
        return (r > 72 ? 72 : r);
    }

    public static int Quadruplex(String str) {
        String Gstruc[] = new String[]{"cccc", "gggg"};
        for (int i = 0; i < Gstruc.length; i++) {
            if (str.contains(Gstruc[i])) {
                return i;
            }
        }
        return -1;
    }

    public static int[] PrimerLocation(String primer, String seq) {
        int n1 = 0;
        int n2 = 0;
        int i = 0;
        int y = 0;
        int x = 0;
        int h = 0;
        int l = 0;
        int z = 0;
        int w = 0;
        int r = 0;

        int k = 7;
        int k1 = 0;
        int k2 = 0;
        int er = 1;
        int sim = 70;
        //  primer = "gcgaaaaccaagtgcttacctcg";

        int l2 = primer.length();
        int l1 = seq.length();
        int lmin = l2;
        if (l1 < 7 || l2 < 7) {
            return null;
        }
        if (lmin > l1) {
            lmin = l1;
        }

        int[] n = new int[0];
        int[] plen = new int[36];
        int[] pmin = new int[36];
        plen[7] = 7;
        pmin[7] = 6;
        plen[8] = 8;
        pmin[8] = 7;
        plen[9] = 9;
        pmin[9] = 7;
        plen[10] = 10;
        pmin[10] = 8;
        plen[11] = 11;
        pmin[11] = 8;
        plen[12] = 12;
        pmin[12] = 9;
        plen[13] = 13;
        pmin[13] = 9;
        plen[14] = 13;
        pmin[14] = 10;
        plen[15] = 13;
        pmin[15] = 10;
        plen[16] = 14;
        pmin[16] = 11;
        plen[17] = 14;
        pmin[17] = 11;
        plen[18] = 15;
        pmin[18] = 12;
        plen[19] = 15;
        pmin[19] = 12;
        plen[20] = 16;
        pmin[20] = 13;
        plen[21] = 16;
        pmin[21] = 13;
        plen[22] = 17;
        pmin[22] = 14;
        plen[23] = 17;
        pmin[23] = 14;
        plen[24] = 18;
        pmin[24] = 15;
        plen[25] = 18;
        pmin[25] = 15;
        plen[26] = 19;
        pmin[26] = 16;
        plen[27] = 19;
        pmin[27] = 16;
        plen[28] = 20;
        pmin[28] = 17;
        plen[29] = 20;
        pmin[29] = 17;
        plen[30] = 21;
        pmin[30] = 18;
        plen[31] = 21;
        pmin[31] = 18;
        plen[32] = 22;
        pmin[32] = 19;
        plen[33] = 22;
        pmin[33] = 19;
        plen[34] = 23;
        pmin[34] = 20;
        plen[35] = 23;
        pmin[35] = 20;
        if (lmin > 6 && lmin < 36) {
            k2 = plen[lmin];
            k1 = pmin[lmin];
        }
        if (lmin > 35) {
            k2 = plen[35];
            k1 = pmin[35];
        }

        int b = 1 + ((k2 * sim) / 100);
        int[] px = new int[l1 + l2];
        char[] o = primer.toCharArray();
        char[] s = seq.toCharArray();

        for (int p = 0; p < l1 - l2 + 1; p++) {
            for (y = 0; y < l2 - k; y++) {
                x = seq.indexOf(primer.substring(l2 - y - k, l2 - y), p);
                if (x == -1) {
                    return n;
                }
                z = x + y + k - 1;
                w = 1 + z - l2;
                l = 0;
                if (px[z] == 0) {
                    n1 = -w;
                    if (n1 < 0) {
                        n1 = 0;
                    }
                    n2 = 1 + z - l1;
                    r = n2;
                    if (n2 < 0) {
                        r = 0;
                        n2 = 0;
                    }
                    if (er - n2 < 0) {
                        px[z] = 1;
                    }
                    if (l2 - n1 - n2 < k2) {
                        px[z] = 1;
                    }
                }
                if (px[z] == 0) {
                    px[z] = 1;
                    for (i = 0; i < k1 - n2; i++) {
                        if ((z - i - n2 < 0) || (l2 - i - n2 - 1 < 0)) {
                            break;
                        }
                        if (s[z - i - n2] == o[l2 - i - n2 - 1]) {
                            l++;
                        } else {
                            r++;
                        }
                    }
                    if (r < er + 1) {
                        for (i = k1 - n2; i < k2 - n2; i++) {
                            if ((z - i - n2 < 0) || (l2 - i - n2 - 1 < 0)) {
                                break;
                            }
                            if (s[z - i - n2] == o[l2 - i - n2 - 1]) {
                                l++;
                            }
                        }
                    }

                    if (l >= b) {
                        h++;
                        if (h == 1) {
                            n = new int[h];
                            n[0] = w + 1;
                        } else {
                            int[] tmp = new int[h];
                            System.arraycopy(n, 0, tmp, 0, h - 1);
                            n = tmp;
                            n[h - 1] = w + 1;
                        }
                        p = w + 1 + l2;
                        break;
                    }
                }
            }
        }
        return n;
    }

    public static int DimersCheck(String primer, String[] listprimer) {
        for (String p : listprimer) {
            if (Oligo.quickDimer(primer, p) > 0) {
                return 1;
            }
        }
        return 0;
    }

}
