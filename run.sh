rm -f nohup.out
CMD='julia --project=. main_mult.jl'
nohup $CMD & disown;
pid=$!
echo julia running in PID=$pid 
ps ux | head -n 1  > sys.log

while true
do
ps ux | grep "$CMD" | grep -v "grep" >> sys.log
sleep 30
done


kill $pid