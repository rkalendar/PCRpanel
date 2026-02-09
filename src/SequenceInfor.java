import java.util.Arrays;

public final class SequenceInfor {

    public SequenceInfor() {
        cdn = new byte[128];
//M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T    
//a   b   c   d   g   h   i   k   m   n   r   s   t   u   v   w   y
//97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 121            
        cdn[65] = 97;   // A
        cdn[66] = 98;   // B
        cdn[67] = 99;   // C
        cdn[68] = 100;  // D
        cdn[71] = 103;  // G
        cdn[72] = 104;  // H
        cdn[73] = 99;   // I
        cdn[75] = 107;  // K
        cdn[77] = 109;  // M
        cdn[78] = 110;  // N
        cdn[82] = 114;  // R
        cdn[83] = 115;  // S
        cdn[84] = 116;  // T
        cdn[85] = 116;  // U
        cdn[86] = 118;  // V
        cdn[87] = 119;  // W
        cdn[89] = 121;  // Y        
        cdn[97] = 97;   // a
        cdn[98] = 98;    // b
        cdn[99] = 99;    // c
        cdn[100] = 100;  // d
        cdn[103] = 103;  // g
        cdn[104] = 104;  // h
        cdn[105] = 99;   // i
        cdn[107] = 107;  // k
        cdn[109] = 109;  // m
        cdn[110] = 110;  // n
        cdn[114] = 114;  // r
        cdn[115] = 115;  // s
        cdn[116] = 116;  // t
        cdn[117] = 116;  // u
        cdn[118] = 118;  // v
        cdn[119] = 119;  // w
        cdn[121] = 121;  // y 
    }

    public String ReadingSequence(String s) {
        byte[] source = s.getBytes();
        if (source == null || source.length == 0) {
            return "";
        }
        dnay = new double[128];
        int n = -1;
        for (int i = 0; i < source.length; i++) {
            if (cdn[source[i]] > 0) {
                dnay[source[i]]++;
                source[++n] = source[i];
            }
        }
        lSeq = n + 1;
        if (n > 1) {
            return new String(Arrays.copyOfRange(source, 0, n + 1));
        }
        return "";
    }

    public double getA() {
        return (dnay[97] + (dnay[109] + dnay[114] + dnay[119]) / 2 + (dnay[118] + dnay[104] + dnay[100]) / 3 + dnay[110] / 4);
    }

    public double getT() {
        return (dnay[116] + dnay[117] + (dnay[121] + dnay[107] + dnay[119]) / 2 + (dnay[98] + dnay[104] + dnay[100]) / 3 + dnay[110] / 4);
    }

    public double getC() {
        return (dnay[99] + (dnay[109] + dnay[115] + dnay[121]) / 2 + (dnay[98] + dnay[118] + dnay[104]) / 3 + dnay[110] / 4);
    }

    public double getG() {
        return (dnay[103] + dnay[105] + ((dnay[115] + dnay[114] + dnay[107]) / 2) + ((dnay[98] + dnay[118] + dnay[100]) / 3) + dnay[110] / 4);
    }

    public double getN() {
        return dnay[110];
    }

    public double getR() {
        return (dnay[103] + dnay[105] + dnay[97] + dnay[114] + ((dnay[115] + dnay[107] + dnay[109] + dnay[119] + dnay[110]) / 2) + ((dnay[98] + dnay[118] + dnay[118] + dnay[104] + dnay[100] + dnay[100]) / 3));
    }

    public double getY() {
        return (dnay[99] + dnay[116] + dnay[117] + dnay[121] + (dnay[109] + dnay[115] + dnay[110] + dnay[107] + dnay[119]) / 2 + (dnay[98] + dnay[118] + dnay[104] + dnay[98] + dnay[104] + dnay[100]) / 3);
    }

    public int getLength() {
        return lSeq;
    }

    public double getCG() {
        return lSeq < 1 ? 0 : (100 * (dnay[103] + dnay[105] + dnay[99] + dnay[115] + ((dnay[114] + dnay[107] + dnay[109] + dnay[121] + dnay[110]) / 2) + ((dnay[98] + dnay[118] + dnay[100] + dnay[98] + dnay[118] + dnay[104]) / 3))) / lSeq;
    }
    private double[] dnay;
    private int lSeq = 0;
    private final byte[] cdn;
}
