import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;

public class OutDirUtil {

    public static Path prepareOutputDir(String outpath) throws IOException {
        if (outpath == null || outpath.isBlank()) {
            throw new IllegalArgumentException("folder_out is empty.");
        }

        Path dir = Paths.get(outpath).toAbsolutePath().normalize();

        // Minimal safeguards against accidental deletion of "the wrong thing"
        guardNotDangerous(dir);

        if (Files.exists(dir)) {
            if (!Files.isDirectory(dir)) {
                throw new IOException("folder_out exists but is not a directory: " + dir);
            }
            clearDirectory(dir);
        } else {
            Files.createDirectories(dir);
        }

        // Check write permission (more reliable than Files.isWritable() on some FS)
        verifyWritable(dir);

        return dir;
    }

    public static void clearDirectory(Path dir) throws IOException {
        try (DirectoryStream<Path> ds = Files.newDirectoryStream(dir)) {
            for (Path child : ds) {
                deleteRecursively(child);
            }
        }
    }

    public static void deleteRecursively(Path path) throws IOException {
        if (!Files.exists(path, LinkOption.NOFOLLOW_LINKS)) return;

        Files.walkFileTree(path, new SimpleFileVisitor<>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                Files.deleteIfExists(file);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                if (exc != null) throw exc;
                Files.deleteIfExists(dir);
                return FileVisitResult.CONTINUE;
            }
        });
    }

    private static void guardNotDangerous(Path dir) throws IOException {
        // dir.getParent()==null оusually means root ("/" on Linux, "C:\" на Windows)
        if (dir.getParent() == null) {
            throw new IOException("Refusing to use a filesystem root as folder_out: " + dir);
        }

        String s = dir.toString();
        if (s.equals(".") || s.equals("..")) {
            throw new IOException("Refusing to use relative special path as folder_out: " + dir);
        }
    }

    private static void verifyWritable(Path dir) throws IOException {
        Path test = dir.resolve(".write_test_" + System.nanoTime() + ".tmp");
        try {
            Files.write(test, "ok".getBytes(StandardCharsets.UTF_8),
                    StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
        } finally {
            Files.deleteIfExists(test);
        }
    }

    public static String stripQuotes(String s) {
        if (s == null) return null;
        s = s.trim();
        if ((s.startsWith("\"") && s.endsWith("\"")) || (s.startsWith("'") && s.endsWith("'"))) {
            if (s.length() >= 2) return s.substring(1, s.length() - 1).trim();
        }
        return s;
    }
}
