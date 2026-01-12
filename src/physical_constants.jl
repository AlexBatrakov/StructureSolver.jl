"""
	G

Newtonian gravitational constant in CGS units.

Unit: cm^3 g^-1 s^-2.
"""
const G =  6.67408e-8


raw"""
	M_sun

Solar mass in CGS units (grams).

Note: computed as $GM_\odot / G$ using `1.3271244e26` (cm^3/s^2) for $GM_\odot$.
"""
const M_sun = 1.3271244e26 / G
#const G = 6.67e-8
"""
	mp

Baryon mass used by the package in CGS units.

Unit: g.
"""
const mp = 1.66e-24
#const mp = 1.66053906660e-24
#const mp = 1.67262192369e-24
"""
	c

Speed of light in vacuum in CGS units.

Unit: cm/s.
"""
const c = 2.99792458e10
const d = 3600*24
const km = 1e5

const k = 1.3800649e-16
const R_gas = 8.31446e7
const σ_st = 5.6704e-5
const a_st = 7.5657e-15