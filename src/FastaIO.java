
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

public class FastaIO {

    /**
     * A simple container for a single FASTA record.
     */
    public static final class FastaRecord {

        public final String id;        // ID  
        public final String header;    // full title (without '>')
        public final String sequence;  // DNA sequence (lowecase, no spaces)

        public FastaRecord(String id, String header, String sequence) {
            this.id = id;
            this.header = header;
            this.sequence = sequence;
        }

        public int length() {
            return sequence.length();
        }

        @Override
        public String toString() {
            return ">" + header + " [length=" + length() + "]";
        }
    }

    /**
     * Stream reading of FASTA records from one or more files. Each file can be
     * plain text or a .gz archive. Calls 'sink.accept(record)' for each parsed
     * record.
     *
     * @param reffiles list of paths to files
     * @param sink handler function for each record

     */
    public static void readFastaStreaming(List<String> reffiles, Consumer<FastaRecord> sink) throws IOException {
        if (reffiles == null || reffiles.isEmpty()) {
            return;
        }

        for (String path : reffiles) {
            try (BufferedReader br = new BufferedReader(new InputStreamReader(openMaybeGzip(path), StandardCharsets.UTF_8), 65536)) { // 64KB буфер

                String line;
                String header = null;
                String id = null;
                StringBuilder sb = new StringBuilder(1048576); // 1 MB 

                int lineNumber = 0;
                while ((line = br.readLine()) != null) {
                    lineNumber++;
                    line = line.trim();

                    if (line.isEmpty()) {
                        continue;
                    }

                    if (line.charAt(0) == '>') {
                        if (header != null) {
                            String seq = dna.DNA(sb.toString());
                            sink.accept(new FastaRecord(id, header, seq));
                            sb.setLength(0);
                        }
                        // Start a new entry
                        header = line.substring(1).trim();
                        id = header.isEmpty() ? "" : header.split("\\s+", 2)[0];
                    } else {
                        if (header == null) {
                            throw new IOException("Error on line " + lineNumber + " в файле " + path + ": sequence found before header '>'");
                        }
                        appendSequence(sb, line);
                    }
                }

                if (header != null) {
                    String seq = dna.DNA(sb.toString());
                    sink.accept(new FastaRecord(id, header, seq));
                }
            }
        }
    }


    /**
     * Convenient method: reads all entries in LinkedHashMap (preserves order).
     * Keys are sequence IDs, values are the sequences themselves.
     *
     * @param reffiles list of paths to files
     * @return map ID -> sequence
     */
    public static Map<String, String> readAllToMap(List<String> reffiles) throws IOException {
        LinkedHashMap<String, String> map = new LinkedHashMap<>();
        readFastaStreaming(reffiles, rec -> {
            String key = rec.id.isEmpty() ? "unnamed" : rec.id;
            // If duplicate IDs are found, add a numerical suffix.
            if (map.containsKey(key)) {
                int i = 2;
                while (map.containsKey(key + "_" + i)) {
                    i++;
                }
                key = key + "_" + i;
            }
            map.put(key, rec.sequence);
        });
        return map;
    }

    /**
     * Reads all entries in the list
     *
     * @param reffiles list of paths to files
     * @return list of all entries
     */
    public static List<FastaRecord> readAllToList(List<String> reffiles) throws IOException {
        List<FastaRecord> records = new ArrayList<>();
        readFastaStreaming(reffiles, records::add);
        return records;
    }

    /**
     * Reads one FASTA record from the file
     *
     * @param filepath path to the file
     * @return first record found or null
     */
    public static FastaRecord readSingle(String filepath) throws IOException {
        final FastaRecord[] holder = new FastaRecord[1];
        readFastaStreaming(Collections.singletonList(filepath), rec -> {
            if (holder[0] == null) {
                holder[0] = rec;
            }
        });
        return holder[0];
    }

    /**
     * Writing FASTA records to a file
     *
     * @param records list of records to write
     * @param outputPath path to the output file
     * @param lineWidth line width for formatting the sequence
     * (usually 60 or 80)
     */
    public static void writeFasta(List<FastaRecord> records, String outputPath, int lineWidth) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(outputPath), StandardCharsets.UTF_8))) {
            for (FastaRecord rec : records) {
                writer.write(">");
                writer.write(rec.header);
                writer.newLine();

                // Write a sequence with line breaks
                String seq = rec.sequence;
                for (int i = 0; i < seq.length(); i += lineWidth) {
                    int end = Math.min(i + lineWidth, seq.length());
                    writer.write(seq, i, end - i);
                    writer.newLine();
                }
            }
        }
    }

    /**
     * Entry with a default line width of 80 characters
     * @param records
     * @param outputPath
     */
    public static void writeFasta(List<FastaRecord> records, String outputPath) throws IOException {
        writeFasta(records, outputPath, 80);
    }

    // === Auxiliary methods ===
    private static InputStream openMaybeGzip(String path) throws IOException {
        InputStream in = new FileInputStream(path);
        if (path.endsWith(".gz")) {
            return new GZIPInputStream(in, 65536); // 64KB буфер для GZIP
        }
        return in;
    }

    private static void appendSequence(StringBuilder sb, String line) {
        for (int i = 0, n = line.length(); i < n; i++) {
            char c = line.charAt(i);
            if (c >= 'A' && c <= 'Z') {
                char u = (char) (c - 32); // toLowerCase for ASCII
                if (u == 'U') {
                    u = 't';
                }
                if (u == 'I') {
                    u = 'g';
                }
                
                sb.append(u);
            } else if (c >= 'a' && c <= 'z') {
                char u = c;
                if (u == 'U') {
                    u = 't';
                }
                                                if (u == 'i') {
                    u = 'g';
                }
                sb.append(u);
            }
        }
    }

    // === Example of use ===
/*
    public static void main(String[] args) {
        try {
            // Пример 1: Чтение одного файла
            System.out.println("=== Чтение FASTA файла ===");
            List<String> files = Arrays.asList("example.fasta");
            
            // Потоковое чтение
            readFastaStreaming(files, record -> {
                System.out.println(record);
                System.out.println("  Последовательность: " + record.sequence.substring(0, Math.min(50, record.sequence.length())) + "...");
                System.out.println();
            });

            // Пример 2: Чтение в список
            System.out.println("\n=== Чтение в список ===");
            List<FastaRecord> records = readAllToList(files);
            System.out.println("Прочитано записей: " + records.size());
            
            // Пример 3: Чтение в карту
            System.out.println("\n=== Чтение в карту ===");
            Map<String, String> seqMap = readAllToMap(files);
            for (String id : seqMap.keySet()) {
                System.out.println(id + " -> длина: " + seqMap.get(id).length());
            }

            // Пример 4: Обратная комплементарная последовательность
            if (!records.isEmpty()) {
                System.out.println("\n=== Reverse Complement ===");
                FastaRecord first = records.get(0);
                String rc = first.reverseComplement();
                System.out.println("Оригинал:  " + first.sequence.substring(0, Math.min(30, first.sequence.length())));
                System.out.println("Rev-Comp:  " + rc.substring(0, Math.min(30, rc.length())));
            }

            // Пример 5: Запись в файл
            System.out.println("\n=== Запись в файл ===");
            writeFasta(records, "output.fasta", 60);
            System.out.println("Записано в output.fasta");

        } catch (IOException e) {
            System.err.println("Ошибка: " + e.getMessage());
            e.printStackTrace();
        }
    }
     */
}
