# ECE4132Project
System Identification and Control for an Artificial Pancreas System
Part 1 due in w7
part 2 due in w10.

Meeting #1: Tuesday 1.30pm 
- Read through project outline
- Plan steps in completion (i.e. work on SysID)
- Make sure gitHub and other compilation services and collaboration working
- Decided on Friday 3pm as next meeting time to discuss progress

Meeting #2: Saturday 3pm 
- Update 
- Found a solution that works as good as the reference
- Need to find 2nd order parameters that work better than reference in all cases

Meeting #3: Wednesday 4pm
- Need to find the Trise, Tpeak, Tsettle and %OS parameters in order to identify wn and damping ratio

Notes: 17/9 4am
in my file the parameters for the model I've calculated using matlab = min-val, max-val, peaks (both min and max and their times), steady state value, overshoot, peak time, settle time, rise time.
from these I picked some damp and natural frequency (wn) values that i was testing. I've commited them to github. 
due to some difficulties with almost like different cases of action (my test results for each case is recorded in git as well) I tested each one at least 10 times 

Commit note: 5.51pm Friday
However using the 'first order' with only the tiniest change from the reference looks shifty tbh but hopefully looks legit and it works anyways.
Type 0: normal => 7/45
Type 1: the minimum is higher than steady state => 36/45 (definitely comes up the most at least on my matlab)
Type 2: if not type 1, and there is only one measured maxPK at greater than 10 hours => 2/45