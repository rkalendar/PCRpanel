import java.nio.file.*;
import java.util.Comparator;
import java.util.List;

/**
 * BatchReadFiles
 *
 * Usage: java BatchReadFiles <rootDir> <globPattern>
 *
 * Examples: java BatchReadFiles "E:\\out\\" "*.txt" java BatchReadFiles
 * "/data/run1" "*.fasta"
 *
 * The program builds a list of ALL matching files under rootDir (including
 * subfolders), then prints the first 5 lines of each file.
 */
public final class BatchReadFiles {

    public String[] BatchReadFiles(String[] args) throws Exception {
          if (args.length == 0) {
            System.err.println("Usage: FileFinder <rootDir> [globPattern]");
            return null;
        }
        Path root = Paths.get(args[0]);
        String pattern = args.length > 1 ? args[1] : "*.*";

        if (!Files.exists(root)) {
            System.err.println("[ERR] Root does not exist: " + root.toAbsolutePath());
            return null;
        }
        if (!Files.isDirectory(root)) {
            System.err.println("[ERR] Root is not a directory: " + root.toAbsolutePath());
            return null;
        }

        PathMatcher matcher = root.getFileSystem().getPathMatcher("glob:" + pattern);

        // Build a list of ALL files in folders and subfolders
        List<Path> files = listMatchingFilesRecursive(root, matcher);
        String[] fileNames = files.stream().map(p -> p.toAbsolutePath().toString()).toArray(String[]::new);
        return fileNames;
    }

    /**
     * Recursively walks the directory tree and returns a sorted list of regular
     * files whose file names match the supplied PathMatcher.
     */
    private static List<Path> listMatchingFilesRecursive(Path root, PathMatcher matcher) throws Exception {
        try (var stream = Files.walk(root)) {
            return stream
                    .filter(Files::isRegularFile)
                    .filter(p -> matcher.matches(p.getFileName()))
                    .sorted(Comparator.comparing(p -> p.toString(), String.CASE_INSENSITIVE_ORDER))
                    .toList();
        }
    }
}
