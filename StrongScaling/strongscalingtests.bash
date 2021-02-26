#!/bin/bash
one=$(mksub strongscaling2.pbs)
echo $one
two=$(mksub -W depend=afterok:$one strongscaling4.pbs)
echo $two
three=$(mksub -W depend=afterok:$two strongscaling8.pbs)
echo $three
four=$(mksub -W depend=afterok:$three strongscaling16.pbs)
echo $four
five=$(mksub -W depend=afterok:$four strongscaling32.pbs)
echo $five
six=$(mksub -W depend=afterok:$five strongscaling64.pbs)
echo $six
seven=$(mksub -W depend=afterok:$six strongscaling128.pbs)
echo $seven
eight=$(mksub -W depend=afterok:$seven strongscaling256.pbs)
echo $eight
