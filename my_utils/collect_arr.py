import pypar

def collect_arr(arr):
    """
    A useful collection routine for pypar.
    If you are using pypar to parallelize the set of nested loops and fill
    the resulting array, you usually need to combine the resulting array
    from several mpi threads. In that case you just can execute
    res=collect_arr(res)
    And it will add the arrays from all the threads and store them in
    the thread number 0
    """
    if pypar.rank() > 0:
        pypar.send(arr, 0)
    else:
        for i in range(1, pypar.size()):
            arr = arr + pypar.receive(i)
    return arr

