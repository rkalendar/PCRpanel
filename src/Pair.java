class Pair {

    int forward;
    int reverse;
    int exonid;
    int pcrsize;

    double rTm;
    double rCG;
    int panelID;         // panel 1-2
    int rln;             // length
    int rlc;             // Linguistic_Complexity
    int rx5;             // location x1
    String rprimer;      // primer sequence
    String rprimername;  // name

    double fTm;
    double fCG;
    int fln;            // length
    int flc;            // Linguistic_Complexity
    int fx5;            // location x1
    String fprimer;     // primer sequence
    String fprimername; // name

    public void SetPanelID(int n) {
        panelID = n;
    }

    public void AddFprimer(String primer, String name, int x5, int ln, int lc, double tm, double cg) {
        this.fprimer = primer;
        this.fprimername = name;
        this.fln = ln;
        this.flc = lc;
        this.fx5 = x5;
        this.fCG = cg;
        this.fTm = tm;
    }

    public void AddRprimer(String primer, String name, int x5, int ln, int lc, double tm, double cg) {
        this.rprimer = primer;
        this.rprimername = name;
        this.rln = ln;
        this.rlc = lc;
        this.rx5 = x5;
        this.rCG = cg;
        this.rTm = tm;
    }

    public Pair(int forward, int reverse, int exonid, int pcrsize) {
        this.forward = forward;
        this.reverse = reverse;
        this.exonid = exonid;
        this.pcrsize = pcrsize;
    }
}
