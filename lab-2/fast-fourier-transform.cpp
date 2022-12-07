#include "./fast-fourier-transform.hpp";


static void drfti1(uint64_t n, double* wa, int* ifac)
{
	static int		ntryh[4] = { 4,2,3,5 };
	static double	tpi = 6.28318530717958647692528676655900577;
	double			arg, argh, argld, fi;
	int64_t				ntry = 0, i, j = -1;
	int64_t				k1, l1, l2, ib;
	int64_t				ld, ii, ip, is, nq, nr;
	int64_t				ido, ipm, nfm1;
	int64_t				nl = n;
	int64_t				nf = 0;

L101:
	++j;

	if (j < 4) ntry = ntryh[j];
	else ntry += 2;

L104:
	nq = nl / ntry;
	nr = nl - ntry * nq;

	if (nr != 0) goto L101;

	nf++;
	ifac[nf + 1] = ntry;
	nl = nq;

	if (ntry != 2) goto L107;
	if (nf == 1) goto L107;

	for (i = 1; i < nf; ++i)
	{
		ib = nf - i + 1;
		ifac[ib + 1] = ifac[ib];
	}
	ifac[2] = 2;

L107:
	if (nl != 1) goto L104;

	ifac[0] = n;
	ifac[1] = nf;

	argh = tpi / n;
	is = 0;
	nfm1 = nf - 1;
	l1 = 1;

	if (nfm1 == 0) return;

	for (k1 = 0; k1 < nfm1; ++k1)
	{
		ip = ifac[k1 + 2];
		ld = 0;
		l2 = l1 * ip;
		ido = n / l2;
		ipm = ip - 1;

		for (j = 0; j < ipm; ++j)
		{
			ld += l1;
			i = is;
			argld = (double)ld * argh;
			fi = 0.0;
			for (ii = 2; ii < ido; ii += 2)
			{
				fi += 1.0;
				arg = fi * argld;
				wa[i++] = cos(arg);
				wa[i++] = sin(arg);
			}
			is += ido;
		}

		l1 = l2;
	}
}

void initFastFourierTransform(uint64_t n, double* wsave, int* ifac)
{
	if (n == 1) return;

	drfti1(n, wsave + n, ifac);
}

static void dradf2(int64_t ido, int64_t l1, double* cc, double* ch, double* wa1)
{
	uint64_t	i, k;
	double		ti2, tr2;
	int64_t		t0, t1, t2, t3, t4, t5, t6;

	t1 = 0;
	t0 = (t2 = l1 * ido);
	t3 = ido << 1;

	for (k = 0; k < l1; ++k) 
	{
		ch[t1 << 1] = cc[t1] + cc[t2];
		ch[(t1 << 1) + t3 - 1] = cc[t1] - cc[t2];
		t1 += ido;
		t2 += ido;
	}

	if (ido < 2) return;
	if (ido == 2) goto L105;

	t1 = 0;
	t2 = t0;

	for (k = 0; k < l1; ++k)
	{
		t3 = t2;
		t4 = (t1 << 1) + (ido << 1);
		t5 = t1;
		t6 = t1 + t1;
		
		for (i = 2; i < ido; i += 2) 
		{
			t3 += 2;
			t4 -= 2;
			t5 += 2;
			t6 += 2;
			tr2 = wa1[i - 2] * cc[t3 - 1] + wa1[i - 1] * cc[t3];
			ti2 = wa1[i - 2] * cc[t3] - wa1[i - 1] * cc[t3 - 1];
			ch[t6] = cc[t5] + ti2;
			ch[t4] = ti2 - cc[t5];
			ch[t6 - 1] = cc[t5 - 1] + tr2;
			ch[t4 - 1] = cc[t5 - 1] - tr2;
		}
		t1 += ido;
		t2 += ido;
	}

	if (ido % 2 == 1) return;

L105:
	t3 = (t2 = (t1 = ido) - 1);
	t2 += t0;
	
	for (k = 0; k < l1; ++k) 
	{
		ch[t1] = -cc[t2];
		ch[t1 - 1] = cc[t3];
		t1 += ido << 1;
		t2 += ido;
		t3 += ido;
	}
}

static void dradf4(int64_t ido, int64_t l1, double* cc, double* ch, double* wa1, double* wa2, double* wa3)
{
	static double	hsqt2 = .70710678118654752440084436210485;
	int64_t			i, k, t0, t1, t2, t3, t4, t5, t6;
	double			ci2, ci3, ci4, cr2, cr3, cr4;
	double			ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

	t0 = l1 * ido;
	t1 = t0;
	t4 = t1 << 1;
	t2 = t1 + (t1 << 1);
	t3 = 0;

	for (k = 0; k < l1; ++k)
	{
		tr1 = cc[t1] + cc[t2];
		tr2 = cc[t3] + cc[t4];
		ch[t5 = t3 << 2] = tr1 + tr2;
		ch[(ido << 2) + t5 - 1] = tr2 - tr1;
		ch[(t5 += (ido << 1)) - 1] = cc[t3] - cc[t4];
		ch[t5] = cc[t2] - cc[t1];

		t1 += ido;
		t2 += ido;
		t3 += ido;
		t4 += ido;
	}

	if (ido < 2) return;
	if (ido == 2) goto L105;

	t1 = 0;
	for (k = 0; k < l1; ++k)
	{
		t2 = t1;
		t4 = t1 << 2;
		t5 = (t6 = ido << 1) + t4;

		for (i = 2; i < ido; i += 2)
		{
			t3 = (t2 += 2);
			t4 += 2;
			t5 -= 2;

			t3 += t0;
			cr2 = wa1[i - 2] * cc[t3 - 1] + wa1[i - 1] * cc[t3];
			ci2 = wa1[i - 2] * cc[t3] - wa1[i - 1] * cc[t3 - 1];
			t3 += t0;
			cr3 = wa2[i - 2] * cc[t3 - 1] + wa2[i - 1] * cc[t3];
			ci3 = wa2[i - 2] * cc[t3] - wa2[i - 1] * cc[t3 - 1];
			t3 += t0;
			cr4 = wa3[i - 2] * cc[t3 - 1] + wa3[i - 1] * cc[t3];
			ci4 = wa3[i - 2] * cc[t3] - wa3[i - 1] * cc[t3 - 1];

			tr1 = cr2 + cr4;
			tr4 = cr4 - cr2;
			ti1 = ci2 + ci4;
			ti4 = ci2 - ci4;
			ti2 = cc[t2] + ci3;
			ti3 = cc[t2] - ci3;
			tr2 = cc[t2 - 1] + cr3;
			tr3 = cc[t2 - 1] - cr3;


			ch[t4 - 1] = tr1 + tr2;
			ch[t4] = ti1 + ti2;

			ch[t5 - 1] = tr3 - ti4;
			ch[t5] = tr4 - ti3;

			ch[t4 + t6 - 1] = ti4 + tr3;
			ch[t4 + t6] = tr4 + ti3;

			ch[t5 + t6 - 1] = tr2 - tr1;
			ch[t5 + t6] = ti1 - ti2;
		}
		t1 += ido;
	}

	if (ido % 2 == 1)return;

L105:

	t2 = (t1 = t0 + ido - 1) + (t0 << 1);
	t3 = ido << 2;
	t4 = ido;
	t5 = ido << 1;
	t6 = ido;

	for (k = 0; k < l1; k++) {
		ti1 = -hsqt2 * (cc[t1] + cc[t2]);
		tr1 = hsqt2 * (cc[t1] - cc[t2]);
		ch[t4 - 1] = tr1 + cc[t6 - 1];
		ch[t4 + t5 - 1] = cc[t6 - 1] - tr1;
		ch[t4] = ti1 - cc[t1 + t0];
		ch[t4 + t5] = ti1 + cc[t1 + t0];
		t1 += ido;
		t2 += ido;
		t4 += t3;
		t6 += ido;
	}
}

static void dradfg(int64_t ido, int64_t ip, int64_t l1, int64_t idl1, double* cc, double* c1, double* c2, double* ch, double* ch2, double* wa)
{
	static double	tpi = 6.28318530717958647692528676655900577;
	int64_t				idij, ipph, i, j, k, l, ic, ik, is;
	int64_t				t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	double			dc2, ai1, ai2, ar1, ar2, ds2;
	int64_t				nbd;
	double			dcp, arg, dsp, ar1h, ar2h;
	int64_t				idp2, ipp2;

	arg = tpi / (double)ip;
	dcp = cos(arg);
	dsp = sin(arg);
	ipph = (ip + 1) >> 1;
	ipp2 = ip;
	idp2 = ido;
	nbd = (ido - 1) >> 1;
	t0 = l1 * ido;
	t10 = ip * ido;

	if (ido == 1) goto L119;

	for (ik = 0; ik < idl1; ++ik)
	{
		ch2[ik] = c2[ik];
	}

	t1 = 0;

	for (j = 1; j < ip; ++j) 
	{
		t1 += t0;
		t2 = t1;
		
		for (k = 0; k < l1; ++k)
		{
			ch[t2] = c1[t2];
			t2 += ido;
		}
	}

	is = -ido;
	t1 = 0;

	if (nbd > l1)
	{
		for (j = 1; j < ip; ++j)
		{
			t1 += t0;
			is += ido;
			t2 = -ido + t1;

			for (k = 0; k < l1; ++k)
			{
				idij = is - 1;
				t2 += ido;
				t3 = t2;
				
				for (i = 2; i < ido; i += 2)
				{
					idij += 2;
					t3 += 2;
					ch[t3 - 1] = wa[idij - 1] * c1[t3 - 1] + wa[idij] * c1[t3];
					ch[t3] = wa[idij - 1] * c1[t3] - wa[idij] * c1[t3 - 1];
				}
			}
		}
	}

	else 
	{
		for (j = 1; j < ip; ++j)
		{
			is += ido;
			idij = is - 1;
			t1 += t0;
			t2 = t1;
			for (i = 2; i < ido; i += 2) {
				idij += 2;
				t2 += 2;
				t3 = t2;
				for (k = 0; k < l1; k++) {
					ch[t3 - 1] = wa[idij - 1] * c1[t3 - 1] + wa[idij] * c1[t3];
					ch[t3] = wa[idij - 1] * c1[t3] - wa[idij] * c1[t3 - 1];
					t3 += ido;
				}
			}
		}
	}

	t1 = 0;
	t2 = ipp2 * t0;

	if (nbd < l1)
	{
		for (j = 1; j < ipph; ++j)
		{
			t1 += t0;
			t2 -= t0;
			t3 = t1;
			t4 = t2;
			
			for (i = 2; i < ido; i += 2)
			{
				t3 += 2;
				t4 += 2;
				t5 = t3 - ido;
				t6 = t4 - ido;

				for (k = 0; k < l1; ++k)
				{
					t5 += ido;
					t6 += ido;
					c1[t5 - 1] = ch[t5 - 1] + ch[t6 - 1];
					c1[t6 - 1] = ch[t5] - ch[t6];
					c1[t5] = ch[t5] + ch[t6];
					c1[t6] = ch[t6 - 1] - ch[t5 - 1];
				}
			}
		}
	}
	else
	{
		for (j = 1; j < ipph; ++j)
		{
			t1 += t0;
			t2 -= t0;
			t3 = t1;
			t4 = t2;
			
			for (k = 0; k < l1; ++k)
			{
				t5 = t3;
				t6 = t4;
				
				for (i = 2; i < ido; i += 2)
				{
					t5 += 2;
					t6 += 2;
					c1[t5 - 1] = ch[t5 - 1] + ch[t6 - 1];
					c1[t6 - 1] = ch[t5] - ch[t6];
					c1[t5] = ch[t5] + ch[t6];
					c1[t6] = ch[t6 - 1] - ch[t5 - 1];
				}

				t3 += ido;
				t4 += ido;
			}
		}
	}

L119:
	for (ik = 0; ik < idl1; ++ik)
	{
		c2[ik] = ch2[ik];
	}

	t1 = 0;
	t2 = ipp2 * idl1;

	for (j = 1; j < ipph; ++j)
	{
		t1 += t0;
		t2 -= t0;
		t3 = t1 - ido;
		t4 = t2 - ido;
		
		for (k = 0; k < l1; ++k)
		{
			t3 += ido;
			t4 += ido;
			c1[t3] = ch[t3] + ch[t4];
			c1[t4] = ch[t4] - ch[t3];
		}
	}

	ar1 = 1.0;
	ai1 = 0.0;
	t1 = 0;
	t2 = ipp2 * idl1;
	t3 = (ip - 1) * idl1;

	for (l = 1; l < ipph; ++l)
	{
		t1 += idl1;
		t2 -= idl1;
		ar1h = dcp * ar1 - dsp * ai1;
		ai1 = dcp * ai1 + dsp * ar1;
		ar1 = ar1h;
		t4 = t1;
		t5 = t2;
		t6 = t3;
		t7 = idl1;

		for (ik = 0; ik < idl1; ++ik)
		{
			ch2[t4++] = c2[ik] + ar1 * c2[t7++];
			ch2[t5++] = ai1 * c2[t6++];
		}

		dc2 = ar1;
		ds2 = ai1;
		ar2 = ar1;
		ai2 = ai1;

		t4 = idl1;
		t5 = (ipp2 - 1) * idl1;
		
		for (j = 2; j < ipph; ++j)
		{
			t4 += idl1;
			t5 -= idl1;

			ar2h = dc2 * ar2 - ds2 * ai2;
			ai2 = dc2 * ai2 + ds2 * ar2;
			ar2 = ar2h;

			t6 = t1;
			t7 = t2;
			t8 = t4;
			t9 = t5;
			
			for (ik = 0; ik < idl1; ++ik)
			{
				ch2[t6++] += ar2 * c2[t8++];
				ch2[t7++] += ai2 * c2[t9++];
			}
		}
	}

	t1 = 0;

	for (j = 1; j < ipph; ++j)
	{
		t1 += idl1;
		t2 = t1;

		for (ik = 0; ik < idl1; ++ik)
		{
			ch2[ik] += c2[t2++];
		}
	}

	if (ido < l1) goto L132;

	t1 = 0;
	t2 = 0;

	for (k = 0; k < l1; ++k)
	{
		t3 = t1;
		t4 = t2;

		for (i = 0; i < ido; i++)
		{
			cc[t4++] = ch[t3++];
		}

		t1 += ido;
		t2 += t10;
	}

	goto L135;

L132:
	for (i = 0; i < ido; ++i) 
	{
		t1 = i;
		t2 = i;

		for (k = 0; k < l1; ++k)
		{
			cc[t2] = ch[t1];
			t1 += ido;
			t2 += t10;
		}
	}

L135:
	t1 = 0;
	t2 = ido << 1;
	t3 = 0;
	t4 = ipp2 * t0;

	for (j = 1; j < ipph; ++j)
	{

		t1 += t2;
		t3 += t0;
		t4 -= t0;

		t5 = t1;
		t6 = t3;
		t7 = t4;

		for (k = 0; k < l1; ++k)
		{
			cc[t5 - 1] = ch[t6];
			cc[t5] = ch[t7];
			t5 += t10;
			t6 += ido;
			t7 += ido;
		}
	}

	if (ido == 1) return;
	if (nbd < l1) goto L141;

	t1 = -ido;
	t3 = 0;
	t4 = 0;
	t5 = ipp2 * t0;

	for (j = 1; j < ipph; ++j) 
	{
		t1 += t2;
		t3 += t2;
		t4 += t0;
		t5 -= t0;
		t6 = t1;
		t7 = t3;
		t8 = t4;
		t9 = t5;

		for (k = 0; k < l1; ++k)
		{
			for (i = 2; i < ido; i += 2)
			{
				ic = idp2 - i;
				cc[i + t7 - 1] = ch[i + t8 - 1] + ch[i + t9 - 1];
				cc[ic + t6 - 1] = ch[i + t8 - 1] - ch[i + t9 - 1];
				cc[i + t7] = ch[i + t8] + ch[i + t9];
				cc[ic + t6] = ch[i + t9] - ch[i + t8];
			}

			t6 += t10;
			t7 += t10;
			t8 += ido;
			t9 += ido;
		}
	}

	return;

L141:

	t1 = -ido;
	t3 = 0;
	t4 = 0;
	t5 = ipp2 * t0;

	for (j = 1; j < ipph; ++j)
	{
		t1 += t2;
		t3 += t2;
		t4 += t0;
		t5 -= t0;

		for (i = 2; i < ido; i += 2)
		{
			t6 = idp2 + t1 - i;
			t7 = i + t3;
			t8 = i + t4;
			t9 = i + t5;

			for (k = 0; k < l1; ++k)
			{
				cc[t7 - 1] = ch[t8 - 1] + ch[t9 - 1];
				cc[t6 - 1] = ch[t8 - 1] - ch[t9 - 1];
				cc[t7] = ch[t8] + ch[t9];
				cc[t6] = ch[t9] - ch[t8];
				t6 += t10;
				t7 += t10;
				t8 += ido;
				t9 += ido;
			}
		}
	}
}

static void drftf1(uint64_t n, double* c, double* ch, double* wa, int* ifac)
{
	uint64_t i, k1, l1, l2;
	int64_t		na, kh, nf;
	int64_t		ip, iw, ido, idl1, ix2, ix3;

	nf = ifac[1];
	na = 1;
	l2 = n;
	iw = n;

	for (k1 = 0; k1 < nf; ++k1)
	{
		kh = nf - k1;
		ip = ifac[kh + 1];
		l1 = l2 / ip;
		ido = n / l2;
		idl1 = ido * l1;
		iw -= (ip - 1) * ido;
		na = 1 - na;

		if (ip != 4) goto L102;

		ix2 = iw + ido;
		ix3 = ix2 + ido;

		if (na != 0) dradf4(ido, l1, ch, c, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1);
		else dradf4(ido, l1, c, ch, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1);

		goto L110;

	L102:
		if (ip != 2) goto L104;
		if (na != 0) goto L103;

		dradf2(ido, l1, c, ch, wa + iw - 1);
		goto L110;

	L103:
		dradf2(ido, l1, ch, c, wa + iw - 1);
		goto L110;

	L104:
		if (ido == 1) na = 1 - na;
		if (na != 0) goto L109;

		dradfg(ido, ip, l1, idl1, c, c, c, ch, ch, wa + iw - 1);
		na = 1;
		goto L110;

	L109:
		dradfg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa + iw - 1);
		na = 0;

	L110:
		l2 = l1;
	}

	if (na == 1)
		return;

	for (i = 0; i < n; ++i)
	{
		c[i] = ch[i];
	}
}

std::shared_ptr<double[]> executeFastFourierTransform(randomSequence* sequence, const uint64_t& bitSize)
{
	int	ifac[15];

	auto bitSequence = std::unique_ptr<uint8_t[]>(sequence->getBitSequence());

	auto X = std::shared_ptr<double[]>(new double[bitSize]);
	auto wsave = std::shared_ptr<double[]>(new double[2 * bitSize]);

	for (uint64_t i = 0; i < bitSize; ++i)
	{
		X[i] = 2 * bitSequence[i] - 1;
	}

	initFastFourierTransform(bitSize, wsave.get(), ifac);
	drftf1(bitSize, X.get(), wsave.get(), wsave.get() + bitSize, ifac);

	return X;
}
