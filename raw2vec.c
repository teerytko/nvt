/*
 * Copyright (c) 2015 Intel Corporation. All rights Reserved.
 */

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <error.h>

#define CLAMP(x,mi,ma)	((x) < (mi) ? (mi) : ((x) > (ma) ? (ma) : (x)))

// BXT A0
static void test_infra_prepare_isys_line_interleaved_image(void *in, int image_width, int image_height, int stride_in, void *out)
{
	static const unsigned int TWO_LINES = 2;
	static const unsigned int TWO_PIXEL_STRIDE = 2;

	int i, j, stride;
	uint32_t *frame_evenline, *frame_oddline;
	uint32_t *in_ptr = (uint32_t *)in;
	uint32_t *out_ptr = (uint32_t *)out;
	stride_in /= 2;						/* Convert from bytes to pixels */

	stride = stride_in/TWO_PIXEL_STRIDE;				/*stride for 32 bit elements */
	frame_oddline = in_ptr;
	frame_evenline = in_ptr + stride;

	for(i=0; i < image_height; i+=2) {
		for(j=0; j < (image_width/TWO_PIXEL_STRIDE); j++) {	/* go over twice the width, 2 lines */
			*out_ptr++ = *(frame_oddline+j);
			*out_ptr++ = *(frame_evenline+j);
		}
		frame_oddline += TWO_LINES * stride;	/* point to odd lines */
		frame_evenline += TWO_LINES * stride;	/* point to even lines */
	}
}

// BXT A0
static void devectorize(void *in, int width, int height, int stride, void *out)
{
	static const unsigned int BYP = 2;	/* Bytes per pixel */
	unsigned int y, x;
	unsigned char *src = in;
	unsigned char *dst = out;

	for (y = 0; y < height; y++) {
		unsigned char *d = dst;
		for (x = 0; x < width; x++) {
			unsigned int ofs = ((x & ~1) << 1) + (x & 1) + ((y & 1) << 1);
			unsigned char *s = &src[y/2*(2*stride) + ofs*BYP];
			d[0] = s[0];
			d[1] = s[1];
			d += BYP;
		}
		dst += stride;
	}
}

// BXT B0
static void test_infra_prepare_bayer_vcc_image(void *in, int image_width, int image_height, int stride_in, void *out)
{
	static const unsigned int BITS_PER_PIXEL = 16;
	static const unsigned int VEC_UINT32S = 32 / 2;

	int row, col;
	int vec_reg_idx  = 0;   // integer to vector indexator
	int frame_padding = 1;
	unsigned int ddr_in_offset = 0;
	unsigned int ddr_out_offset = 0;
	uint16_t bayer_c00_01[4] = {0};
	uint16_t bayer_c10_11[4] = {0};
	uint32_t vec_buf_c00[VEC_UINT32S]; // Bayer vectors composition space
	uint32_t vec_buf_c01[VEC_UINT32S];
	uint32_t vec_buf_c10[VEC_UINT32S];
	uint32_t vec_buf_c11[VEC_UINT32S];
	uint32_t p0,p1;
	unsigned char *p_in = (unsigned char *)in;
	unsigned char *p_out = (unsigned char *)out;
	uint32_t vec_stride = 64;

	stride_in /= 2;						/* Convert from bytes to pixels */
	memset(vec_buf_c00, 0, sizeof(vec_buf_c00));
	memset(vec_buf_c01, 0, sizeof(vec_buf_c01));
	memset(vec_buf_c10, 0, sizeof(vec_buf_c10));
	memset(vec_buf_c11, 0, sizeof(vec_buf_c11));

	// some assumptions
	if (image_width % 2 != 0)
		error(1, 0, "ERROR: incorrect image width\n");
	if (image_height %2 != 0)
		error(1, 0, "ERROR: incorrect image height\n");

	// determine stride (need multiple of 64)
	if (image_width % 64) {
		// need padding
		printf("#### WARNING: not yet focussed on padding, double check\n");
	}
	// now fill the buffer
	// Bayer vector building - send 32 components of each color component
	ddr_out_offset = 0;
	for (row = 0 ; row < image_height ; row += 2) {
		for (col = 0; col < image_width ; col += 4) {
			// load components from DDR
			ddr_in_offset = (row * stride_in * sizeof(uint16_t)) + (col * sizeof(uint16_t));
			//buffer_load(in + ddr_in_offset, &bayer_c00_01, sizeof(bayer_c00_01));
			memcpy(&bayer_c00_01, p_in + ddr_in_offset, sizeof(bayer_c00_01));
			ddr_in_offset = ((row + 1) * stride_in * sizeof(uint16_t)) + (col * sizeof(uint16_t));
			//buffer_load(in + ddr_in_offset, &bayer_c10_11, sizeof(bayer_c10_11));
			memcpy(&bayer_c10_11, p_in + ddr_in_offset, sizeof(bayer_c10_11));
			// process vectors
			p0 = bayer_c00_01[0];
			p1 = bayer_c00_01[2];
			vec_buf_c00[vec_reg_idx] = (p0 | (p1 << BITS_PER_PIXEL));
			p0 = bayer_c00_01[1];
			p1 = bayer_c00_01[3];
			vec_buf_c01[vec_reg_idx] = (p0 | (p1 << BITS_PER_PIXEL));
			p0 = bayer_c10_11[0];
			p1 = bayer_c10_11[2];
			vec_buf_c10[vec_reg_idx] = (p0 | (p1 << BITS_PER_PIXEL));
			p0 = bayer_c10_11[1];
			p1 = bayer_c10_11[3];
			vec_buf_c11[vec_reg_idx] = (p0 | (p1 << BITS_PER_PIXEL));
			vec_reg_idx++;
			if (0 == (vec_reg_idx % VEC_UINT32S))  // The 32 pixels of each plain are ready to be sent
			{
				//buffer_store(out + ddr_out_offset, &vec_buf_c00, sizeof(vec_buf_c00));
				memcpy(p_out + ddr_out_offset, &vec_buf_c00, sizeof(vec_buf_c00));
				ddr_out_offset += vec_stride;
				//buffer_store(out + ddr_out_offset, &vec_buf_c01, sizeof(vec_buf_c01));
				memcpy(p_out + ddr_out_offset, &vec_buf_c01, sizeof(vec_buf_c01));
				ddr_out_offset += vec_stride;
				//buffer_store(out + ddr_out_offset, &vec_buf_c10, sizeof(vec_buf_c10));
				memcpy(p_out + ddr_out_offset, &vec_buf_c10, sizeof(vec_buf_c10));
				ddr_out_offset += vec_stride;
				//buffer_store(out + ddr_out_offset, &vec_buf_c11, sizeof(vec_buf_c11));
				memcpy(p_out + ddr_out_offset, &vec_buf_c11, sizeof(vec_buf_c11));
				ddr_out_offset += vec_stride;
				// Re-use the buffer (all the up-till-now built vectors have been sent)
				vec_reg_idx = 0;
			}
		}
	}
	// write out last set of vectors if input image is not big enough
	if (vec_reg_idx && frame_padding) {
		//buffer_store(out + ddr_out_offset, &vec_buf_c00, sizeof(vec_buf_c00));
		memcpy(p_out + ddr_out_offset, &vec_buf_c00, sizeof(vec_buf_c00));
		ddr_out_offset += vec_stride;
		//buffer_store(out + ddr_out_offset, &vec_buf_c01, sizeof(vec_buf_c01));
		memcpy(p_out + ddr_out_offset, &vec_buf_c01, sizeof(vec_buf_c01));
		ddr_out_offset += vec_stride;
		//buffer_store(out + ddr_out_offset, &vec_buf_c10, sizeof(vec_buf_c10));
		memcpy(p_out + ddr_out_offset, &vec_buf_c10, sizeof(vec_buf_c10));
		ddr_out_offset += vec_stride;
		//buffer_store(out + ddr_out_offset, &vec_buf_c11, sizeof(vec_buf_c11));
		memcpy(p_out + ddr_out_offset, &vec_buf_c11, sizeof(vec_buf_c11));
		ddr_out_offset += vec_stride;
		// Re-use the buffer (all the up-till-now built vectors have been sent)
		vec_reg_idx = 0;
	}
}

// BXT B0
static void devectorize_b0(void *in, int width, int height, int stride, void *out)
{
	static const unsigned int BYP = 2;	/* Bytes per pixel */
	unsigned int y, x, p;
	unsigned char *src = in;
	unsigned char *dst = out;

	if (width & 1) error(1, 0, "ERROR: incorrect image width\n");
	if (height & 1) error(1, 0, "ERROR: incorrect image height\n");

	for (y = 0; y < height; y += 2) {
		unsigned char *d = dst;
		unsigned char *s = src;
		for (x = 0; x < width; x += 64) {
//if (x==5312/*-64*/)
			for (p = 0; p < 32; p++) {
				if (x + p*2 >= width) break;
				d[p*4+0] = s[p*2+0];			// Gr
				d[p*4+1] = s[p*2+1];
				d[p*4+2] = s[p*2+64];			// R
				d[p*4+3] = s[p*2+65];
				d[stride+p*4+0] = s[p*2+128];		// B
				d[stride+p*4+1] = s[p*2+129];
				d[stride+p*4+2] = s[p*2+192];		// Gb
				d[stride+p*4+3] = s[p*2+193];
			}
			s += 4 * 32 * 2;
			d += 64 * BYP;
		}
		dst += 2 * stride;
		src += 2 * stride;
	}
}

// Input:  Tiled-Y
// Output: V4L2_PIX_FMT_NV12
static void devectorize_tiley(void *in, int width, int height, int stride, void *out)
{
	static const unsigned int BLOCK_W = 16;
	static const unsigned int BLOCK_H = 32;
	const unsigned int BLOCK_SIZE = BLOCK_W * BLOCK_H;
	const unsigned int block_cols = (width + BLOCK_W - 1) / BLOCK_W;
	const unsigned int block_rows = (height + BLOCK_H - 1) / BLOCK_H;
	// Luma plane addresses for source and destination
	unsigned char *s = in;
	unsigned char *d = out;
	// Chroma plane addresses for source and destination
	unsigned char *sc = &s[BLOCK_SIZE * block_cols * block_rows];
	unsigned char *dc = &d[height * stride];
	int x, y;

	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			int block_x = x / BLOCK_W;
			int block_y = y / BLOCK_H;
			int block_ofs_x = x % BLOCK_W;
			int block_ofs_y = y % BLOCK_H;

			unsigned char *block = &s[BLOCK_SIZE * (block_cols * block_y + block_x)];
			d[y * stride + x] = block[block_ofs_y * BLOCK_W + block_ofs_x];

			if ((x & 1) == 0 && (y & 1) == 0) {
				// Handle chroma
				int chroma_x = x / 2;
				int chroma_y = y / 2;

				block_y = chroma_y / BLOCK_H;
				block_ofs_x = chroma_x % (BLOCK_W / 2);
				block_ofs_y = chroma_y % BLOCK_H;
				block = &sc[BLOCK_SIZE * (block_cols * block_y + block_x)];
				dc[chroma_y * stride + chroma_x * 2]     = block[block_ofs_y * BLOCK_W + block_ofs_x * 2];
				dc[chroma_y * stride + chroma_x * 2 + 1] = block[block_ofs_y * BLOCK_W + block_ofs_x * 2 + 1];
			}
		}
	}
}

// Input:  V4L2_PIX_FMT_YYUV420_V32
// Output: V4L2_PIX_FMT_NV12
static void devectorize_yuv(void *in, int width, int height, int stride, void *out)
{
	static const int VEC_SIZE = (64 + 2*32 + 64) * 2;	// in bytes
	static const int LUMA_SHIFT = 8;	// In theory 8, 7 gives brighter image
	static const int CHROMA_SHIFT = 6;	// Should be verified
	const int w = width;
	const int h = height;
	int p, y, x, y0, x0, c;
	unsigned char *d = out;
	unsigned char *s = in;

	if (width & 63) printf("WARNING: width would be better to be multiple of 64");
	if (height & 1) error(1, 0, "height must be multiple of 2");

	for (y = 0; y < h; y+= 2) {
		for (x = 0; x < w; x += 64) {
			// luma
			for (y0 = 0; y0 < 2; y0++) for (x0 = 0; x0 < 64; x0++) {
				uint16_t *s0 = (uint16_t *)s;
				int x1 = (x0 & 31) | (y0 << 5);
				int y1 = (x0 & 32) >> 5;
				p = s0[y1 * 128 + x1] >> LUMA_SHIFT;
				d[(y+y0)*stride + x+x0] = CLAMP(p, 0, 255);
			}
			// chroma
			for (x0 = 0; x0 < 32; x0++) for (c = 0; c < 2; c++) {
				int16_t *s0 = (int16_t *)s;
				p = s0[64 + x0 + 32 * c];
				d[stride*h + stride*(y/2) + x + x0*2 + c] =
						CLAMP((p >> CHROMA_SHIFT) + 128, 0, 255);
			}
			s += VEC_SIZE;
		}
	}
}

static void *read_file(const char *filename, int *size)
{
	int r;
	void *p;
	FILE *f;

	f = fopen(filename, "rb");              if (!f) error(1, 0, "failed to open file");
	r = fseek(f, 0, SEEK_END);              if (r < 0) error(1, 0, "failed to seek file");
	*size = ftell(f);                       if (*size < 0) error(1, 0, "failed to get file size");
	r = fseek(f, 0, SEEK_SET);              if (r < 0) error(1, 0, "failed to seek file");
	p = calloc(1, *size);                   if (!p) error(1, 0, "out of memory");
	r = fread(p, *size, 1, f);              if (r != 1) error(1, 0, "failed to read file");
	fclose(f);

	return p;
}

static void write_file(const char *name, const unsigned char *data, int size)
{
	FILE *f;
	int r;

	f = fopen(name, "wb");			if (!f)	error(1, 0, "can not open file");
	r = fwrite(data, size, 1, f);		if (r != 1) error(1, 0, "failed to write data to file");
	r = fclose(f);				if (r != 0) error(1, 0, "failed to close file");
}

int main(int argc, char *argv[])
{
	int width, height, stride;
	int size, out_size;
	void *in, *out;
	int b0;
	int tiley;
	int vecy;
	int devec;

	if (argc != 6) {
		printf("Usage: %s <bayvd> <width> <height> <in file> <out file>\n", argv[0]);
		printf( "First argument:\n"
			"b - vectorized raw Bayer BXT B0 format\n"
			"a - vectorized raw Bayer BXT A0 format\n"
			"y - vectorized tile-y format\n"
			"v - vectorized yuv format (V4L2_PIX_FMT_YYUV420_V32)\n"
			"d - devectorize from the format (default vectorize)\n"
			"Output is either the vectorized format or (with d)\n"
			"devectorized raw (bayer or NV24) format\n");
		exit(1);
	}

	b0 = strchr(argv[1], 'b') != NULL;
	tiley = strchr(argv[1], 'y') != NULL;
	vecy = strchr(argv[1], 'v') != NULL;
	devec = strchr(argv[1], 'd') != NULL;

	width = atoi(argv[2]);
	height = atoi(argv[3]);
	in = read_file(argv[4], &size);

	if (tiley) {
		stride = width;
		if (size != stride*height*3/2) { printf("WARNING: read %i bytes, expected %i\n", size, stride*height*3/2); }
		out_size = width * height * 3;
		out = calloc(1, out_size);

		if (!devec) {
			error(1, 0, "conversion to tile-y not supported");
		} else {
			devectorize_tiley(in, width, height, stride, out);
		}
	} else if (vecy) {
		int width64 = (width + 63) & ~63;
		stride = width;
		if (size != width64 * height * 3) { printf("WARNING: read %i bytes, expected %i\n", size, width64 * height * 3); }
		out_size = width * height * 3 / 2;
		out = calloc(1, out_size);

		if (!devec) {
			error(1, 0, "conversion to vector-y not supported");
		} else {
			devectorize_yuv(in, width, height, stride, out);
		}
	} else {
		printf("TEEMUR: in final else\n");
		stride = width * 2;
		if (size != stride*height) { printf("WARNING: read %i bytes, expected %i\n", size, stride*height); }
		out_size = width * height * 2;
		out = calloc(1, out_size);

		if (!b0) {
			if (!devec) {
				test_infra_prepare_isys_line_interleaved_image(in, width, height, stride, out);
			} else {
				devectorize(in, width, height, stride, out);
			}
		} else {
			if (!devec) {
				test_infra_prepare_bayer_vcc_image(in, width, height, stride, out);
			} else {
				devectorize_b0(in, width, height, stride, out);
			}
		}
	}
	write_file(argv[5], out, out_size);

	return 0;
}
