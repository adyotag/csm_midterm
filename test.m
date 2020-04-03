% Benchmark for time

a = Assembly("Biaxial_Q4_2x2.txt").run();
b = Assembly("Beam_Bending_Q4_4x1_Al.txt").run();
c = Assembly("Beam_Bending_Q8_4x1_Al.txt").run();
d = Assembly("Beam_Bending_Q9_4x1_Al.txt").run();
e = Assembly("Beam_Bending_Q4_8x2_Al.txt").run();
f = Assembly("Beam_Bending_Q4_16x4_Al.txt").run();
g = Assembly("Beam_Bending_Q4_16x4_PU.txt").run();
h = Assembly("Beam_Bending_Q9_16x4_PU.txt").run();