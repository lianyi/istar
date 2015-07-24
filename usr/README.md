To compile from source, one needs to change the OPENBABEL_ROOT variable in the Makefile, to point to the directory where openbabel can be found, e.g. ~/openbabel-2.3.2

The calculations of the three moments are as folllows:

	m1 = sum(di) / n

	m2 = sqrt(sum((di - m1) ^ 2) / n)

	m3 = cbrt(sum((di - m1) ^ 3) / n)

All the three moments are of unit Armstrong.

These three formulae are also implemented in both [USR:OptIso] by Ting Zhou et al. and [ElectroShape] by M. Stuart Armstrong et al.

Note that the calculations of m2 and m3 in the first version of [USR] do not compute sqrt() and cbrt().

Also note that the calculation of m3 in [USRCAT] is slightly different:

	m3 = cbrt(sum((di - m1) ^ 3) / n) / m2 = cbrt(skew(di))

This calculation has three disadvantages:

1) Overflow could occur when m2 = 0, e.g. in the case that there is only one atom in one of the five subsets. In this case, USRCAT uses a branch statement to output three zeros for (m1, m2, m3).

2) m3 is unitless, while m1 and m2 are of unit Armstrong.

3) computationally slower because of an extra division.

[USR:OptIso]: http://dx.doi.org/10.1016/j.jmgm.2010.08.007
[ElectroShape]: http://dx.doi.org/10.1007/s10822-011-9463-8
[USR]: http://dx.doi.org/10.1002/jcc.20681
[USRCAT]: http://dx.doi.org/10.1186/1758-2946-4-27
