# Instructions:
In circuit.txt file put
(total number of nodes)
Then list every element in the circuit in the following format
(node from) (node to) (element type) (element value)
For voltage controlled sources
(node from) (node to) (element type) (multiplier) (controlNodeFrom) (controlNodeTo)

##(element type):
0 for uncontrolled current source (node from is current entering source, node to is current exiting source)
1 for uncontrolled voltage source (node from is negative, node to is positive)
2 for resistor (node from is negative terminal, node to is positive terminal)
3 for voltage controlled voltage source
4 for voltage controlled current source

eg.
4 --> 4 nodes
1 2 1 20 --> from node 1 to node 2, uncontrolled voltage source of 20 Volts (negative @ node 1, positive @ node 2)
3 2 0 5 --> from node 3 to node 2, uncontrolled current source of 5 Amps (current going from node 3 to node 2)
3 1 2 4 --> from node 3 to node 1, resistor of 4 Ohms
4 2 3 0.8 1 3 --> from node 4 to node 2, voltage controlled voltage source with voltage 0.8(V3-V1)
2 3 4 2.1 3 4 --> from node 2 to node 3, current controlled voltage source with current 2.1(V4-V3)

Don't need to pick a reference node or pre-calculate nodes connected to voltage sources (the program will do it for you)
Note that there will be a free parameter because the voltage of the reference node is arbitrary
Therefore, solution determined by program may have all node voltages offset by a constant compared to 
    solution by hand if a different reference node is picked.
Take caution not to repeat any elements (they will be treated as 2 identical elements in parallel)
Cannot solve circuits with current controlled sources
Solved circuit will be outputted to results.txt
Might add mesh current analysis soon
Have fun!