
# ./longrun_ve.sh | tee ./longrun_ve_out.txt

rm -rf test3
rm -rf test3_le_?
rm -rf test3_ve_?

#RUNPATH='/anaconda3/envs/pde20191119/bin/'                    # mac
#RUNPATH='/home/icsrsss/anaconda3/envs/blockdata/bin/'         # server (11,12)
RUNPATH='/usr/bin/'                                            # docker

date
echo -e  '20\n 10\n   4000\n 1' | $RUNPATH/python simulator.py
mv -f test3 test3_ve_1
echo; echo ' - - - - - - E N D   O  F   R U N  - - - - - - - - - - - - - -'; echo
date
echo -e  '40\n 20\n   8000\n 1' | $RUNPATH/python simulator.py 
mv -f test3 test3_ve_2
echo; echo ' - - - - - - E N D   O  F   R U N  - - - - - - - - - - - - - -'; echo
date
echo -e  '60\n 30\n  12000\n 1' | $RUNPATH/python simulator.py 
mv -f test3 test3_ve_3
echo; echo ' - - - - - - E N D   O  F   R U N  - - - - - - - - - - - - - -'; echo
date
echo -e  '80\n 40\n  16000\n 1' | $RUNPATH/python simulator.py 
mv -f test3 test3_ve_4
echo; echo ' - - - - - - E N D   O  F   R U N  - - - - - - - - - - - - - -'; echo
date
echo -e '100\n 50\n  20000\n 1' | $RUNPATH/python simulator.py 
mv -f test3 test3_ve_5
echo; echo ' - - - - - - E N D   O  F   R U N  - - - - - - - - - - - - - -'; echo
date
echo -e '120\n 60\n  24000\n 1' | $RUNPATH/python simulator.py 
mv -f test3 test3_ve_6
echo; echo ' - - - - - - E N D   O  F   R U N  - - - - - - - - - - - - - -'; echo


date
