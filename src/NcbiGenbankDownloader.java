import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.URLEncoder;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.time.Duration;
import java.util.*;
import java.util.stream.Collectors;

public class NcbiGenbankDownloader {

    private static final String EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

    private final String tool;
    private final String email;
    private final String apiKey; // можно null
    private final HttpClient http;

    public NcbiGenbankDownloader(String tool, String email, String apiKey) {
        this.tool = Objects.requireNonNull(tool);
        this.email = Objects.requireNonNull(email);
        this.apiKey = apiKey;
        this.http = HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(20)).build();
    }

    public static void main(String[] args) throws Exception {
        String geneSymbol = "BRCA2";
        String taxId = "9606";              // Homo sapiens
        Integer ngFrom = 13732;             // example range (can be null)
        Integer ngTo = 58896;               //  example range (can be null)

        var dl = new NcbiGenbankDownloader(
                "my_java_ncbi_tool",
                "youremail@example.org",
                System.getenv("NCBI_API_KEY") // optional
        );

        Path outDir = Paths.get("out", geneSymbol);
        Files.createDirectories(outDir);

        Optional<String> geneId = dl.findGeneId(geneSymbol, taxId);
        if (geneId.isEmpty()) {
            System.err.println("GeneID not found for " + geneSymbol);
            return;
        }
        System.out.println("GeneID=" + geneId.get());

        // --- NM_ only (RefSeq RNA) ---
        List<String> refseqRnaAccs = dl.elinkAccessionVersions(geneId.get(), "gene_nuccore_refseqrna");
        List<String> nmAccs = refseqRnaAccs.stream()
                .filter(a -> a.startsWith("NM_"))
                .distinct()
                .toList();

        System.out.println("NM_ records: " + nmAccs.size());
        for (String acc : nmAccs) {
            Path out = outDir.resolve(acc + ".gb");
            dl.efetchGenbank(acc, out, null, null); // полный GenBank
            dl.politeDelay();
        }

        // --- RefSeqGene genomic (NG_) ---
        List<String> refseqGeneAccs = dl.elinkAccessionVersions(geneId.get(), "gene_nuccore_refseqgene");
        List<String> ngAccs = refseqGeneAccs.stream()
                .filter(a -> a.startsWith("NG_"))
                .distinct()
                .toList();

        System.out.println("NG_ (RefSeqGene) records: " + ngAccs.size());
        for (String acc : ngAccs) {
            Path out = outDir.resolve(acc + ".gb");
            dl.efetchGenbank(acc, out, ngFrom, ngTo); // range as from/to (1-based), or null/null for full
            dl.politeDelay();
        }
    }

    // ---------- Step 1: GeneID ----------
    public Optional<String> findGeneId(String geneSymbol, String taxId) throws Exception {
        // символ + organism (txid)
        String term = geneSymbol + "[Gene Name] AND txid" + taxId + "[Organism]";
        URI uri = uri("esearch.fcgi", Map.of(
                "db", "gene",
                "term", term,
                "retmode", "xml",
                "retmax", "5"
        ));
        Document doc = getXml(uri);
        List<String> ids = xpathText(doc, "//IdList/Id/text()");
        return ids.isEmpty() ? Optional.empty() : Optional.of(ids.get(0));
    }

    // ---------- Step 2: ELink gene->nuccore, сразу accession.version ----------
    public List<String> elinkAccessionVersions(String geneId, String linkname) throws Exception {
        // idtype=acc => ELink returns accession.version :contentReference[oaicite:6]{index=6}
        URI uri = uri("elink.fcgi", Map.of(
                "dbfrom", "gene",
                "db", "nuccore",
                "id", geneId,
                "linkname", linkname,
                "idtype", "acc",
                "retmode", "xml"
        ));
        Document doc = getXml(uri);
        // В ELink XML это будут <Link><Id>ACC.V</Id></Link>
        return xpathText(doc, "//LinkSetDb/Link/Id/text()");
    }

    // ---------- Step 3: EFetch GenBank (gbwithparts) + (optional) range ----------
    public void efetchGenbank(String accver, Path outFile, Integer seqStart, Integer seqStop) throws Exception {
        // GenBank flat file: rettype=gb или gbwithparts; we take gbwithparts :contentReference[oaicite:7]{index=7}
        Map<String, String> p = new LinkedHashMap<>();
        p.put("db", "nuccore");
        p.put("id", accver);
        p.put("rettype", "gbwithparts");
        p.put("retmode", "text");

        // диапазон 1-based (как "from/to"): seq_start/seq_stop :contentReference[oaicite:8]{index=8}
        if (seqStart != null && seqStop != null) {
            p.put("seq_start", String.valueOf(seqStart));
            p.put("seq_stop", String.valueOf(seqStop));
        }

        URI uri = uri("efetch.fcgi", p);

        HttpRequest req = HttpRequest.newBuilder(uri)
                .timeout(Duration.ofMinutes(2))
                .GET()
                .build();

        HttpResponse<InputStream> resp = http.send(req, HttpResponse.BodyHandlers.ofInputStream());
        if (resp.statusCode() != 200) throw new IOException("EFetch HTTP " + resp.statusCode() + " for " + accver);

        try (InputStream in = resp.body();
             OutputStream out = Files.newOutputStream(outFile, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)) {
            in.transferTo(out);
        }
    }

    // ---------- helpers ----------
    private void politeDelay() {
        try { Thread.sleep((apiKey == null || apiKey.isBlank()) ? 400 : 150); } catch (InterruptedException ignored) {}
    }

    private URI uri(String endpoint, Map<String, String> params) {
        Map<String, String> p = new LinkedHashMap<>(params);
        // NCBI usage policy: tool + email; api_key при высокой частоте :contentReference[oaicite:9]{index=9}
        p.put("tool", tool);
        p.put("email", email);
        if (apiKey != null && !apiKey.isBlank()) p.put("api_key", apiKey);

        String qs = p.entrySet().stream()
                .map(e -> enc(e.getKey()) + "=" + enc(e.getValue()))
                .collect(Collectors.joining("&"));
        return URI.create(EUTILS + endpoint + "?" + qs);
    }

    private static String enc(String s) {
        return URLEncoder.encode(s, StandardCharsets.UTF_8);
    }

    private Document getXml(URI uri) throws Exception {
        HttpRequest req = HttpRequest.newBuilder(uri).timeout(Duration.ofSeconds(30)).GET().build();
        HttpResponse<InputStream> resp = http.send(req, HttpResponse.BodyHandlers.ofInputStream());
        if (resp.statusCode() != 200) throw new IOException("HTTP " + resp.statusCode() + " for " + uri);

        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        dbf.setFeature("http://apache.org/xml/features/disallow-doctype-decl", true);
        dbf.setExpandEntityReferences(false);

        try (InputStream in = resp.body()) {
            return dbf.newDocumentBuilder().parse(in);
        }
    }

    private static List<String> xpathText(Document doc, String expr) throws Exception {
        var xp = XPathFactory.newInstance().newXPath();
        NodeList nl = (NodeList) xp.evaluate(expr, doc, XPathConstants.NODESET);
        List<String> out = new ArrayList<>(nl.getLength());
        for (int i = 0; i < nl.getLength(); i++) out.add(nl.item(i).getNodeValue());
        return out;
    }
}
