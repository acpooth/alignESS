# Log de pruebas

## Pruebas

| Comando | Resultado |
|--------|:-----------:|
| alignESS.py pair 2.7.1:5.3.1:5.3.1:2.7.1:4.1.2:1.2.1 5.3.1:5.3.1:4.2.1 | OK |
| alignESS.py -l  pair 2.7.1:5.3.1:5.3.1:2.7.1:4.1.2:1.2.1 5.3.1:5.3.1:4.2.1 | OK |
| alignESS.py dbalign  nr\_part.db | OK |
| alignESS.py dbalign -db2 random\_frag.txt   nr\_part.db | OK | 
| alignESS.py dbalign -db2 multi.txt   nr\_part.db | [Error][1] |
| alignESS.py dbalign -db2 random_frag.txt -nproc 4 -align -t 0.7   nr\_part.db | OK |
| alignESS.py multi multi.txt | OK | 
| test\_nwx.py | [Warning][2] 

## Errores 

[1]: Error relacionado a que __multi.txt__ tiene texto antes de los códigos nùmericos 

``` Bash
------ Opening databases:
nr_part.db ( sqlite )
multi.txt ( text_noind )
------ Aligning both databases ------
Number of processes: 2
multiprocessing.pool.RemoteTraceback: 
"""
Traceback (most recent call last):
  File "/home/cperalta/.conda/envs/ess-env/lib/python3.10/multiprocessing/pool.py", line 125, in worker
    result = (True, func(*args, **kwds))
  File "/home/cperalta/.conda/envs/ess-env/lib/python3.10/multiprocessing/pool.py", line 48, in mapstar
    return list(map(*args))
  File "nw_ec_alignx.pyx", line 476, in nw_ec_alignx.seq_vs_db
  File "nw_ec_alignx.pyx", line 367, in nw_ec_alignx.NW
  File "nw_ec_alignx.pyx", line 90, in nw_ec_alignx._FastNW
KeyError: '#Cluster14'
"""

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/home/cperalta/Documents/repos/alignESS/tests/../alignESS.py", line 480, in <module>
    main()
  File "/home/cperalta/Documents/repos/alignESS/tests/../alignESS.py", line 469, in main
    main_db(args)
  File "/home/cperalta/Documents/repos/alignESS/tests/../alignESS.py", line 390, in main_db
    scores = nwx.db_vs_db(db1, db2, hmat, decs, thres=args.threshold,
  File "nw_ec_alignx.pyx", line 517, in nw_ec_alignx.db_vs_db
    resd = pool.map(align_func, seqs1)
  File "/home/cperalta/.conda/envs/ess-env/lib/python3.10/multiprocessing/pool.py", line 364, in map
    return self._map_async(func, iterable, mapstar, chunksize).get()
  File "/home/cperalta/.conda/envs/ess-env/lib/python3.10/multiprocessing/pool.py", line 771, in get
    raise self._value
KeyError: '#Cluster14'
```
[2]: No output 
