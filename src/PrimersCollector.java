public class PrimersCollector {

    private double[] Tm;
    private double[] CG;
    private int[] ln; // length
    private int[] pl; // Linguistic_Complexity
    private int[] x1;  // location x1

    private int taskid; // task number
    private String[] primer; // primer sequence
    private String[] primername; //name
    private int x1t; // task x1 coordinate
    private int x2t;       // task end coordinate
    private int n;         // number of entrance

    public int Amount() {
        return n;
    }

    public String[] getPrimer() {
        return primer;
    }

    public String[] getpPrimerName() {
        return primername;
    }

    public int[] getPrimerLengths() {
        return pl;
    }

    public int[] getPrimerLC() {
        return ln;
    }

    public int[] getPrimerLocations() {
        return x1;
    }

    public int getTaskId() {
        return taskid;
    }

    public void setTaskId(int t) {
        taskid = t;
    }

    public void setTaskX1(int x1) {
        x1t = x1;
    }

    public void setTaskX2(int x1) {
        x2t = x1;
    }

    public int getTaskX1() {
        return x1t;
    }

    public int getTaskX2() {
        return x2t;
    }

    public double[] getTms() {
        return Tm;
    }

    public double[] getCGs() {
        return CG;
    }

    public void insert(int n, String[] primer, String[] primername, int[] pl, double[] Tm, double[] cg, int[] ln, int[] p) {
        this.n = n;
        this.primer = primer;
        this.primername = primername;
        this.Tm = Tm;
        this.CG = cg;
        this.ln = ln;
        this.pl = pl;
        this.x1 = p;
    }
}
