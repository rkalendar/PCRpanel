import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;

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

public class NcbiGenBankFetcher {

    // NCBI E-utilities base
    private static final String EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

    // Please provide your identifiers (NCBI recommends tool+email)
    private final String tool;
    private final String email;
    private final String apiKey; // optional (can null)

    private final HttpClient http;

    public NcbiGenBankFetcher(String tool, String email, String apiKey) {
        this.tool = Objects.requireNonNull(tool);
        this.email = Objects.requireNonNull(email);
        this.apiKey = apiKey;
        this.http = HttpClient.newBuilder()
                .connectTimeout(Duration.ofSeconds(20))
                .build();
    }

    public static void main(String[] args) throws Exception {
        //  Homo sapiens (TaxID 9606)
        String taxId = "9606";
        List<String> genes = List.of("TP53", "BRCA1", "EGFR", "KRAS");

        // substitute your values:
        String tool = "my_java_fetcher";
        String email = "youremail@example.org";
        String apiKey = System.getenv("NCBI_API_KEY"); // can be omitted

        NcbiGenBankFetcher f = new NcbiGenBankFetcher(tool, email, apiKey);

        Path outDir = Paths.get("out");
        Files.createDirectories(outDir);

        for (String gene : genes) {
            System.out.println("== " + gene + " ==");

            Optional<String> geneId = f.findGeneId(gene, taxId);
            if (geneId.isEmpty()) {
                System.out.println("  GeneID not found");
                continue;
            }
            System.out.println("  GeneID: " + geneId.get());

            // 1) RefSeq RNA (usually mRNA/transcripts)
            List<String> nuccoreIds = f.linkGeneToNuccore(geneId.get(), "gene_nuccore_refseqrna");

            // fallback: RefSeqGene (genomic "reference gene" region)
            if (nuccoreIds.isEmpty()) {
                nuccoreIds = f.linkGeneToNuccore(geneId.get(), "gene_nuccore_refseqgene");
            }

            // fallback: any nuccore
            if (nuccoreIds.isEmpty()) {
                nuccoreIds = f.linkGeneToNuccore(geneId.get(), "gene_nuccore");
            }

            if (nuccoreIds.isEmpty()) {
                System.out.println("  No linked nuccore IDs");
                continue;
            }
            System.out.println("  nuccore IDs: " + nuccoreIds.size());

            Path outFile = outDir.resolve(gene + ".gb");
            f.fetchGenBankNuccore(nuccoreIds, outFile, "gbwithparts"); // или "gb"
            System.out.println("  saved: " + outFile.toAbsolutePath());

            // We observe the limits (without a key ~3 rps; with a key, it is possible to do it more often).
            f.politeDelay();
        }
    }

    /** Search for GeneID by symbol and TaxID. */
    public Optional<String> findGeneId(String geneSymbol, String taxId) throws Exception {
        // Field tags for Gene include "Gene Name" (symbol/aliases) :contentReference[oaicite:4]{index=4}
        // For organisms, we use txidNNNN[Organism] — the most stable option
        String term = geneSymbol + "[Gene Name] AND txid" + taxId + "[Organism]";

        URI uri = uri("esearch.fcgi", Map.of(
                "db", "gene",
                "term", term,
                "retmode", "xml",
                "retmax", "5"
        ));

        Document doc = getXml(uri);
        List<String> ids = x(doc, "//IdList/Id/text()");
        return ids.isEmpty() ? Optional.empty() : Optional.of(ids.get(0));
    }

    /** Link GeneID → nuccore (with linkname specified).
     * @param geneId
     * @param linkname
     * @return 
     * @throws java.lang.Exception */
    public List<String> linkGeneToNuccore(String geneId, String linkname) throws Exception {
        URI uri = uri("elink.fcgi", Map.of(
                "dbfrom", "gene",
                "db", "nuccore",
                "id", geneId,
                "linkname", linkname,
                "retmode", "xml"
        ));
        Document doc = getXml(uri);
        // We take the Id from LinkSetDb (you can filter more precisely by LinkName)
        return x(doc, "//LinkSetDb/Link/Id/text()");
    }

    /** EFetch GenBank flatfile for nuccore id-shnikov. */
    public void fetchGenBankNuccore(List<String> nuccoreIds, Path outFile, String rettype) throws Exception {
        // EFetch rettype=gb or gbwithparts and retmode=text :contentReference[oaicite:5]{index=5}
        // Important: id can be passed through a comma
        String idParam = String.join(",", nuccoreIds);

        URI uri = uri("efetch.fcgi", Map.of(
                "db", "nuccore",
                "id", idParam,
                "rettype", rettype,
                "retmode", "text"
        ));

        HttpRequest req = HttpRequest.newBuilder(uri)
                .GET()
                .timeout(Duration.ofMinutes(2))
                .build();

        HttpResponse<InputStream> resp = http.send(req, HttpResponse.BodyHandlers.ofInputStream());
        if (resp.statusCode() != 200) {
            throw new IOException("EFetch failed HTTP " + resp.statusCode());
        }

        try (InputStream in = resp.body();
             OutputStream out = Files.newOutputStream(outFile, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)) {
            in.transferTo(out);
        }
    }

    /** A slight delay to avoid exceeding the limits. */
    private void politeDelay() {
        try {
// without apiKey, ~350-400 ms between requests is better (≈2.5–3 rps)
// with apiKey, it can be less, but it's also better to be careful
            Thread.sleep(apiKey == null || apiKey.isBlank() ? 400 : 150);
        } catch (InterruptedException ignored) {}
    }

    // ----------------------- helpers -----------------------

    private URI uri(String endpoint, Map<String, String> params) {
        Map<String, String> p = new LinkedHashMap<>(params);
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
        HttpRequest req = HttpRequest.newBuilder(uri)
                .GET()
                .timeout(Duration.ofSeconds(30))
                .build();

        HttpResponse<InputStream> resp = http.send(req, HttpResponse.BodyHandlers.ofInputStream());
        if (resp.statusCode() != 200) {
            throw new IOException("HTTP " + resp.statusCode() + " for " + uri);
        }

        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        // XML security settings
        dbf.setFeature("http://apache.org/xml/features/disallow-doctype-decl", true);
        dbf.setExpandEntityReferences(false);

        try (InputStream in = resp.body()) {
            return dbf.newDocumentBuilder().parse(in);
        }
    }

    private static List<String> x(Document doc, String xpathExpr) throws Exception {
        XPath xp = XPathFactory.newInstance().newXPath();
        NodeList nl = (NodeList) xp.evaluate(xpathExpr, doc, XPathConstants.NODESET);
        List<String> out = new ArrayList<>();
        for (int i = 0; i < nl.getLength(); i++) out.add(nl.item(i).getNodeValue());
        return out;
    }
}
