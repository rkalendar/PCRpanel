public final class DNA2 {

    private int lseq;
    private int[] cdn = new int[128];
    private double[] r = new double[128];

    public DNA2(String seq) {
        if (seq == null || seq.isEmpty()) {
            return;
        }
        cdn = new int[128];
        r = new double[128];
        
        cdn[65] = 2;  // A
        cdn[66] = 2;  // B
        cdn[67] = 2;  // C
        cdn[68] = 2;  // D
        cdn[71] = 2;  // G
        cdn[72] = 2;  // H
        cdn[73] = 2;  // I
        cdn[75] = 2;  // K
        cdn[77] = 2;  // M
        cdn[78] = 2;  // N
        cdn[82] = 2;  // R
        cdn[83] = 2;  // S
        cdn[84] = 2;  // T
        cdn[85] = 2;  // U
        cdn[86] = 2;  // V
        cdn[87] = 2;  // W
        cdn[89] = 2;  // Y        
        cdn[97] = 1;   // a
        cdn[98] = 1;   // b
        cdn[99] = 1;   // c
        cdn[100] = 1;  // d
        cdn[103] = 1;  // g
        cdn[104] = 1;  // h
        cdn[105] = 1;  // i
        cdn[107] = 1;  // k
        cdn[109] = 1;  // m
        cdn[110] = 1;  // n
        cdn[114] = 1;  // r
        cdn[115] = 1;  // s
        cdn[116] = 1;  // t
        cdn[117] = 1;  // u
        cdn[118] = 1;  // v
        cdn[119] = 1;  // w
        cdn[121] = 1;  // y  
        
        int l = seq.length();
        byte[] b = seq.getBytes();
        int t = 0;
        lseq = 0;
        for (int i = 0; i < l; i++) {
            if (b[i] < 122 && b[i] > 64) {
                t = cdn[b[i]];
                if (t > 0) {
                    lseq++;
                    if (t == 1) {
                        r[b[i]]++;
                    } else {
                        r[b[i] + 32]++;
                    }
                }
            }
        }

    }
    public double getCG() {
        return (lseq < 1) ? 0 : ((r[103] + r[99] + r[115] + r[105] + (r[114] + r[109] + r[121] + r[107] + r[110]) / 2 + (r[98] + r[98] + r[118] + r[118] + r[100] + r[104]) / 3) * 100 / lseq);
    }
    public double getTm65() {
        return (lseq < 1) ? 0 : (64.9 + (41 * (r[103] + r[99] + r[115] + r[105] + (r[114] + r[109] + r[121] + r[107] + r[110]) / 2 + (r[98] + r[98] + r[118] + r[118] + r[100] + r[104]) / 3 - 16.4)) / lseq);
    }

}
