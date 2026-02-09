import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

public class VariantGenerator {

    private int stln = 0;
    private final String[] stg;

    public VariantGenerator(String pol) {
        stg = GeneratorVariants(pol);
    }

    public String generatedstr() {
        return String.join(" ", stg);
    }

    public int minstrlen() {
        return stln;
    }

    private String[] GeneratorVariants(String pol) {
        if (pol == null) {
            return null;
        }
        if ("".equals(pol)) {
            String[] d = new String[1];
            d[0] = "";
            stln = 0;
            return d;
        }

        List<String> b = new ArrayList<>();

        String[] parts = pol.split("/");
        for (String part : parts) {
            int k = 1;
            pol = dna.DNA(part);
            stln = pol.length();
            for (int i = 0; i < pol.length(); i++) {
                String s = strgen(pol.charAt(i));
                String[] c = s.split(" ");
                int t = c.length;
                if (t > 1) {
                    k = k * t;
                }
            }
            String[] d = new String[k];
            Arrays.fill(d, "");

            int f = 1;
            for (int i = 0; i < pol.length(); i++) {
                String s = strgen(pol.charAt(i));
                String[] c = s.split(" ");
                int t = c.length;
                if (t > 1) {
                    List<String> z = new ArrayList<>();
                    for (int j = 0; j < f; j++) {
                        for (int m = 0; m < t; m++) {
                            z.add(d[j] + c[m]);
                        }
                    }
                    f = z.size();
                    for (int j = 0; j < f; j++) {
                        d[j] = z.get(j);
                    }
                    z.clear();
                } else {
                    for (int j = 0; j < f; j++) {
                        d[j] += pol.charAt(i);
                    }
                }
            }
            b.addAll(Arrays.asList(d));
        }
        return b.toArray(String[]::new);
    }

    private String strgen(char input) {
        Map<Character, String> mappings = new HashMap<>();
        mappings.put('b', "t c g");
        mappings.put('d', "a t g");
        mappings.put('h', "a t c");
        mappings.put('k', "t g");
        mappings.put('m', "a c");
        mappings.put('n', "a t c g");
        mappings.put('r', "a g");
        mappings.put('s', "c g");
        mappings.put('v', "a c g");
        mappings.put('w', "a t");
        mappings.put('y', "t c");

        return mappings.getOrDefault(input, String.valueOf(input));
    }

}
