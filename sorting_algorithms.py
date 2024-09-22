# Bubble Sort
def bubble_sort(L):
    """
    Sorts a list using the Bubble Sort algorithm.

    Parameters:
    L (list): The list to be sorted.

    Returns:
    None: The list is sorted in place.
    """
    length = len(L)
    for i in range(length - 1):
        for j in range(length - i - 1):
            if L[j] > L[j + 1]:
                L[j], L[j + 1] = L[j + 1], L[j]
    # Time complexity: O(n²)


# Insertion Sort
def insertion_sort(L):
    """
    Sorts a list using the Insertion Sort algorithm.

    Parameters:
    L (list): The list to be sorted.

    Returns:
    None: The list is sorted in place.
    """
    for i in range(len(L)):
        j = i
        while j > 0 and L[j] < L[j - 1]:
            L[j], L[j - 1] = L[j - 1], L[j]
            j = j - 1
    # Time complexity: Worst case: O(n²), Best case: O(n)


# Selection Sort
def selection_sort(L):
    """
    Sorts a list using the Selection Sort algorithm.

    Parameters:
    L (list): The list to be sorted.

    Returns:
    None: The list is sorted in place.
    """
    length = len(L)
    for i in range(length):
        for j in range(i + 1, length):
            if L[i] > L[j]:
                L[j], L[i] = L[i], L[j]
    # Time complexity: O(n²)


# Quick Sort
def quick_sort(L):
    """
    Sorts a list using the Quick Sort algorithm.

    Parameters:
    L (list): The list to be sorted.

    Returns:
    list: A sorted version of the input list.
    """
    if not L:
        return L
    pivot = L[0]
    sup = []
    eq = []
    inf = []

    for i in L:
        if i < pivot:
            inf.append(i)
        elif i > pivot:
            sup.append(i)
        else:
            eq.append(i)

    sup = quick_sort(sup)
    inf = quick_sort(inf)

    return inf + eq + sup
    # Time complexity: Best case: O(n log n), Worst case: O(n²)


# Merge Sort
def fusion(L1, L2):
    """
    Merges two sorted lists into one sorted list.

    Parameters:
    L1 (list): The first sorted list.
    L2 (list): The second sorted list.

    Returns:
    list: A merged sorted list.
    """
    n = len(L1)
    m = len(L2)

    if n == 0:
        return L2
    if m == 0:
        return L1

    i, j = 0, 0
    L = []

    while i < n and j < m:
        if L1[i] < L2[j]:
            L.append(L1[i])
            i += 1
        else:
            L.append(L2[j])
            j += 1

    return L + L2[j:] if i == n else L + L1[i:]


def merge_sort(L):
    """
    Sorts a list using the Merge Sort algorithm.

    Parameters:
    L (list): The list to be sorted.

    Returns:
    list: A sorted version of the input list.
    """
    if len(L) <= 1:
        return L

    n = len(L) // 2
    L1 = L[:n]
    L2 = L[n:]

    L1 = merge_sort(L1)
    L2 = merge_sort(L2)

    return fusion(L1, L2)
    # Time complexity: O(n log n)


# Heap Sort
def comparer(x, y):
    """
    Compares two lists element-wise.

    Parameters:
    x (list): The first list.
    y (list): The second list.

    Returns:
    int: 1 if x > y, 0 if x < y.
    """
    for i in range(len(x)):
        if x[i] > y[i]:
            return 1
        elif x[i] < y[i]:
            return 0


def correct_heap(L, deb, end):
    """
    Maintains the heap property in a list.

    Parameters:
    L (list): The list to be heapified.
    deb (int): The index of the current root element.
    end (int): The last index to consider in the heap.

    Returns:
    None: The heap property is maintained in place.
    """
    end = min(end, len(L) - 1)
    p = deb

    while p < end:
        f = 2 * p + 1
        if f > end:
            return
        if f + 1 <= end and comparer(L[f + 1], L[f]) == 1:
            f += 1
        if comparer(L[f], L[p]) == 0:
            return
        L[p], L[f] = L[f], L[p]
        p = f


def make_heap(L):
    """
    Builds a max heap from an unsorted list.

    Parameters:
    L (list): The list to be converted into a heap.

    Returns:
    None: The heap is built in place.
    """
    for p in range((len(L) - 2) // 2, -1, -1):
        correct_heap(L, p, len(L) - 1)


def heap_sort(L):
    """
    Sorts a list using the Heap Sort algorithm.

    Parameters:
    L (list): The list to be sorted.

    Returns:
    None: The list is sorted in place.
    """
    make_heap(L)
    for f in range(len(L) - 1, -1, -1):
        L[0], L[f] = L[f], L[0]
        correct_heap(L, 0, f - 1)
    # Time complexity: O(n log n)


# Radix Sort
def radix_sort(arr):
    """
    Sorts a list of non-negative integers using the Radix Sort algorithm.

    Parameters:
    arr (list): The list of integers to be sorted.

    Returns:
    list: A sorted version of the input list.
    """
    max_num = max(arr)  # Find the maximum number to determine the number of digits
    exp = 1  # Exponent representing the current digit position

    while max_num // exp > 0:
        bins = {i: [] for i in range(10)}  # Create bins for each digit (0 to 9)

        for number in arr:
            digit = (number // exp) % 10
            bins[digit].append(number)

        arr = []
        for i in range(10):
            arr.extend(bins[i])

        exp *= 10

    return arr
    # Time complexity: O(n * k), where k is the number of digits
