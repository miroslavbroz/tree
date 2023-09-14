
# Tree implementation in Fortran or Python

A tree is a structure like this:

```
              [tree%root]                                 
            /              \                              
      [left]                [right]         ---> x, sorted
     /      \              /       \                      
[left]       [right] [left]         [right] ---> y, sorted

```
It is used to find a nearest neighbor like this:

```
  p +  [node.p]        
      /        \       
[left]          [right]
                       
            <------>   
              diff     
                       
---------------------->
          axis = x or y

if p in on the left from [node.p], open [left]
only if diff < best, open [right]!
otherwise, it can hardly be closer
because left-node-right are closest, sorted in x (or y)
```
See test.f90 for a comparison of this nearest-neighbor search and a brute-force.

