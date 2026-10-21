/**
 * Hash table with open addressing.
 * Key: long, Value: int (positive)
 * No boxing, no objects 5-7 times faster than standard HashMap.
 */
public class LongIntHashMap {

    public static final int  NO_VALUE    = -1;
    private static final long EMPTY_KEY   = Long.MIN_VALUE;
    private static final long DELETED_KEY = Long.MIN_VALUE + 1;
    private static final int  MAX_CAPACITY = 1 << 30; // 1_073_741_824

    private int pruneThreshold = 1; // start by deleting singles (value == 1)

    private long[] keys;
    private int[]  values;
    private int    size;
    private int    capacity;
    private int    threshold;

    // ---------- Constructors ----------

    public LongIntHashMap() {
        this(1 << 20); // ~1 M slots, no initial resize needed
    }

    public LongIntHashMap(int initialCapacity) {
        capacity  = nextPowerOfTwo(initialCapacity);
        keys      = new long[capacity];
        values    = new int[capacity];
        threshold = (int) (capacity * 0.75);
        java.util.Arrays.fill(keys, EMPTY_KEY);
    }

    // ---------- Public API ----------

    /**
     * Add / update value. value must be >= 0.
     */
    public void put(long key, int value) {     
        if (size >= threshold) {
            if (capacity < MAX_CAPACITY) {
                resize();
            } else {
                int removed = pruneBelowOrEqual(pruneThreshold);
                if (removed == 0) {
                    pruneThreshold++;
                    removed = pruneBelowOrEqual(pruneThreshold);
                } else {
                    pruneThreshold = 1; // reset after successful cleanup
                }
                if (removed > 0) {
                    rehash();
                }
                if (size >= threshold) {
                    throw new IllegalStateException(
                        "Cannot free space. size=" + size
                        + ", pruneThreshold=" + pruneThreshold);
                }
            }
        }

        int i = index(key);
        while (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY && keys[i] != key) {
            i = next(i);
        }
        if (keys[i] != key) {
            size++;
        }
        keys[i]   = key;
        values[i] = value;
    }

    /**
     * Returns the value for the key, or NO_VALUE (-1) if absent.
     */
    public int get(long key) {
        int i = index(key);
        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                return values[i];
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    /** Returns true if the key exists in the map. */
    public boolean containsKey(long key) {
        return get(key) != NO_VALUE;
    }

    /**
     * Removes the key. Returns the deleted value, or NO_VALUE if absent.
     */
    public int remove(long key) {
        int i = index(key);
        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                int old = values[i];
                keys[i]   = DELETED_KEY; // tombstone
                values[i] = 0;
                size--;
                return old;
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    public int     size()    { return size; }
    public boolean isEmpty() { return size == 0; }

    public void clear() {
        java.util.Arrays.fill(keys, EMPTY_KEY);
        java.util.Arrays.fill(values, 0);
        size           = 0;
        pruneThreshold = 1;
    }

    /**
     * Merges all entries of this map into target, summing values and
     * capping at maxValue. Used to merge per-thread maps.
     */
    public void mergeInto(LongIntHashMap target, int maxValue) {
        for (int i = 0; i < capacity; i++) {
            long k = keys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;
            int existing = target.get(k);
            if (existing == NO_VALUE) {
                target.put(k, Math.min(values[i], maxValue));
            } else {
                target.put(k, Math.min(existing + values[i], maxValue));
            }
        }
    }

    // ---------- Private helpers ----------

    private int index(long key) {
        int h = (int) (key ^ (key >>> 32));
        h ^= (h >>> 16);
        return h & (capacity - 1);
    }

    private int next(int i) {
        return (i + 1) & (capacity - 1);
    }

    /**
     * Doubles the table. 
     */
    private void resize() {
        int    newCapacity = capacity << 1;
        long[] newKeys     = new long[newCapacity];
        int[]  newValues   = new int[newCapacity];
        java.util.Arrays.fill(newKeys, EMPTY_KEY);

        for (int i = 0; i < capacity; i++) {
            long k = keys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;

            int j = (int) (k ^ (k >>> 32));
            j ^= (j >>> 16);
            j &= (newCapacity - 1);
            while (newKeys[j] != EMPTY_KEY) {
                j = (j + 1) & (newCapacity - 1);
            }
            newKeys[j]   = k;
            newValues[j] = values[i];
        }

        capacity  = newCapacity;
        threshold = (int) (capacity * 0.75);
        keys      = newKeys;
        values    = newValues;
    }

    /**
     * Rehashes the current table in-place (same capacity).
     * Called after pruneBelowOrEqual to eliminate tombstone chains
     * and restore O(1) average probe length.
     */
    private void rehash() {
        long[] newKeys   = new long[capacity];
        int[]  newValues = new int[capacity];
        java.util.Arrays.fill(newKeys, EMPTY_KEY);

        for (int i = 0; i < capacity; i++) {
            long k = keys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;

            int j = index(k);
            while (newKeys[j] != EMPTY_KEY) {
                j = next(j);
            }
            newKeys[j]   = k;
            newValues[j] = values[i];
        }

        keys   = newKeys;
        values = newValues;
    }



    /**
     * Removes all entries with value <= minValue.
     * Returns the number of removed entries.
     */
    private int pruneBelowOrEqual(int minValue) {
        int removed = 0;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY
                    && values[i] <= minValue) {
                keys[i]   = DELETED_KEY;
                values[i] = 0;
                size--;
                removed++;
            }
        }
        return removed;
    }

    private static int nextPowerOfTwo(int n) {
        if (n <= 1) return 2;
        if (n > MAX_CAPACITY) {
            throw new IllegalArgumentException(
                "initialCapacity exceeds maximum: " + MAX_CAPACITY);
        }
        int p = Integer.highestOneBit(n);
        return (p == n) ? n : p << 1;
    }

    // ---------- toString ----------

    @Override
    public String toString() {
        StringBuilder sb    = new StringBuilder("{");
        boolean       first = true;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY) {
                if (!first) sb.append(", ");
                sb.append(keys[i]).append('=').append(values[i]);
                first = false;
            }
        }
        return sb.append('}').toString();
    }

    // ---------- Quick test ----------
    /*
        LongIntHashMap map = new LongIntHashMap();

        map.put(1_000_000_000_000_000_001L, 100);
        map.put(9_223_372_036_854_775_000L, 200);
        map.put(42L, 999);

        System.out.println(map.get(42L));                        // 999
        System.out.println(map.get(1_000_000_000_000_000_001L)); // 100
        System.out.println(map.get(999L));                       // -1 (not found)

        System.out.println(map.containsKey(42L));  // true
        System.out.println(map.containsKey(999L)); // false

        map.remove(42L);
        System.out.println(map.containsKey(42L)); // false
        System.out.println(map.size());           // 2
        System.out.println(map);                  // {key=value, ...}
    */
}
