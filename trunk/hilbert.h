#pragma once
//#define HILBERT
inline uint64 i64log2(register uint64 x)
{
	register uint64 l=0UL;
	if(x >= 1ULL<<32) { x>>=32; l|=32; }
	if(x >= 1ULL<<16) { x>>=16; l|=16; }
	if(x >= 1ULL<<8) { x>>=8; l|=8; }
	if(x >= 1ULL<<4) { x>>=4; l|=4; }
	if(x >= 1ULL<<2) { x>>=2; l|=2; }
	if(x >= 1ULL<<1) l|=1;
	return l;
}


inline uint64 encode_hilbert(const uint64 n, uint64 x, uint64 y)
{
	// the resolution less 1 than can hold the numbers x, y
	const auto res = i64log2(max(x,y));
	// swap when parities of n and res+1 are different(or n and res are of the same parity)
	//	if(!((res+n)&1)) //when the sum is even
	if( !( (res&1) ^ (n&1) ))
		swap(x,y);
	auto code = 0ULL;
	register auto w = 1ULL<<res;
	while(w > 0) {
		const auto a = x < w, b = y < w;
		const auto quadrant = 2*int(!a) + int(a^b);
		code = (code<<2) + quadrant;
		switch(quadrant) {
			case 0:
				swap(x,y);
				break;
			case 1:
				y-= w;
				break;
			case 2:
				x-= w;
				y-= w;
				break;
			case 3:
				auto xold = x;
				x = w - y - 1;
				y = 2*w - xold - 1;
				break;
		}
		w>>=1; //w/=2
	}
	return code;
}

inline void decode_hilbert(const uint64 n, register uint64 code, uint64 &xout, uint64 &yout)
{
	// the resolution less 1 than can hold the numbers x, y
	const auto res = (i64log2(code) >> 1) + 1;
	auto quadrant = code & 3;
	uint64 x, y;
	switch(quadrant) {
			case 0:
				x=0;y=0;
				break;
			case 1:
				x=0;y=1;
				break;
			case 2:
				x=1;y=1;
				break;
			case 3:
				x=1;y=0;
				break;
	}
	register auto w = 2ULL;
	code >>= 2;
	while(code>0) {
		const auto quadrant = code & 3;
		switch(quadrant) {
			case 0:
				swap(x,y);
				break;
			case 1:
				y+= w;
				break;
			case 2:
				x+= w;
				y+= w;
				break;
			case 3:
				auto xold = x;
				x = 2*w - y - 1;
				y = w - xold - 1;
				break;
		}
		w<<=1; //w*=2
		code >>= 2;
	}
	// swap when parities of n and res are different
	if( (res&1) ^ (n&1) ) {
		xout = y;
		yout = x;
	} else {
		xout = x;
		yout = y;
	}
}

class Hilbert2DSpaceFillingCurve {
	const uint8 resolution;
public:
	Hilbert2DSpaceFillingCurve(uint64 nx, uint64 ny) :
	  resolution(i64log2(max(nx,ny)) + 1)
	{	}
	const uint64 length() const {return 1ull << (resolution << 1);}
	uint64 encode(uint64 x, uint64 y) {return encode_hilbert(resolution, x, y);}
	void decode(uint64 h, uint64 &xout, uint64 &yout) {decode_hilbert(resolution, h,xout,yout);}
};