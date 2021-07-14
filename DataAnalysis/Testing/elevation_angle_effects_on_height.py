import math

elv_base = 10  # baseline

elv_10 = 12
elv_14 = 15


gate = 50

Re = 6370  # Radius of the Earth, [km]

range_here = 90 + 15 / 2 + 15 * gate  # Slant range [km]

h_base = math.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here * math.sin(math.radians(elv_base))) - Re
h10 = math.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here * math.sin(math.radians(elv_10))) - Re
h14 = math.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here * math.sin(math.radians(elv_14))) - Re

print(h_base)
print(h10)
print(h14)
