#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>

#include "linux/videodev2.h"
#include "linux/atomisp.h"

char *name = "v4l2n";

#define FALSE		0
#define TRUE		(!FALSE)

#define STRINGIFY_(x)	#x
#define STRINGIFY(x)	STRINGIFY_(x)

#define CLEAR(x)	memset(&(x), 0, sizeof(x));
#define SIZE(x)		(sizeof(x)/sizeof((x)[0]))
#define BIT(x)		(1<<(x))

typedef unsigned char bool;

static struct {
	int verbosity;
	int fd;
	struct v4l2_requestbuffers reqbufs;
} vars = {
	.verbosity = 2,
	.fd = -1,
	.reqbufs.count = 2,
	.reqbufs.type = V4L2_BUF_TYPE_VIDEO_CAPTURE,
	.reqbufs.memory = V4L2_MEMORY_USERPTR,
};

struct symbol_list {
	int id;
	const char *symbol;
};
#define SYMBOL_END	{ -1, NULL }

struct token_list {
	char id;
	int flags;
	const char *token;
	const struct symbol_list *symbols;
};
#define TOKEN_END	{ 0, 0, NULL, NULL }
#define TOKEN_F_OPTARG	BIT(0)		/* Token may have optional argument */
#define TOKEN_F_ARG	BIT(1)		/* Token must have argument or it is an error */
#define TOKEN_F_ARG2	BIT(2)		/* Token may have 2 integer arguments */
#define TOKEN_F_ARG4	BIT(3)		/* Token may have 4 integer arguments */

#define V4L2_CID	"V4L2_CID_"
#define CONTROL(id)	{ V4L2_CID_##id, (#id) }
static const struct symbol_list controls[] = {
	CONTROL(BRIGHTNESS),

	/* Flash controls */
	CONTROL(FLASH_DURATION),
	CONTROL(FLASH_INTENSITY),
	CONTROL(TORCH_INTENSITY),
	CONTROL(INDICATOR_INTENSITY),
	CONTROL(FLASH_TRIGGER),
	CONTROL(FLASH_MODE),
	SYMBOL_END
};

#define V4L2_BUF_TYPE	"V4L2_BUF_TYPE_"
#define BUFTYPE(id)	{ V4L2_BUF_TYPE_##id, (#id) }
static const struct symbol_list buf_types[] = {
	BUFTYPE(VIDEO_CAPTURE),
	// BUFTYPE(VIDEO_CAPTURE_MPLANE),
	BUFTYPE(VIDEO_OUTPUT),
	// BUFTYPE(VIDEO_OUTPUT_MPLANE),
	BUFTYPE(VIDEO_OVERLAY),
	BUFTYPE(VBI_CAPTURE),
	BUFTYPE(VBI_OUTPUT),
	BUFTYPE(SLICED_VBI_CAPTURE),
	BUFTYPE(SLICED_VBI_OUTPUT),
	BUFTYPE(VIDEO_OUTPUT_OVERLAY),
	BUFTYPE(PRIVATE),
	SYMBOL_END
};

#define V4L2_PIX_FMT	"V4L2_PIX_FMT_"
#define PIXFMT(id)	{ V4L2_PIX_FMT_##id, (#id) }
static const struct symbol_list pixelformats[] = {
	PIXFMT(RGB332),
	PIXFMT(RGB444),
	PIXFMT(RGB555),
	PIXFMT(RGB565),
	PIXFMT(RGB555X),
	PIXFMT(RGB565X),
	PIXFMT(BGR24),
	PIXFMT(RGB24),
	PIXFMT(BGR32),
	PIXFMT(RGB32),
	PIXFMT(GREY),
	PIXFMT(Y4),
	PIXFMT(Y6),
	PIXFMT(Y10),
	PIXFMT(Y16),
	PIXFMT(PAL8),
	PIXFMT(YVU410),
	PIXFMT(YVU420),
	PIXFMT(YUYV),
	PIXFMT(YYUV),
	PIXFMT(YVYU),
	PIXFMT(UYVY),
	PIXFMT(VYUY),
	PIXFMT(YUV422P),
	PIXFMT(YUV411P),
	PIXFMT(Y41P),
	PIXFMT(YUV444),
	PIXFMT(YUV555),
	PIXFMT(YUV565),
	PIXFMT(YUV32),
	PIXFMT(YUV410),
	PIXFMT(YUV420),
	PIXFMT(HI240),
	PIXFMT(HM12),
	PIXFMT(NV12),
	PIXFMT(NV21),
	PIXFMT(NV16),
	PIXFMT(NV61),
	PIXFMT(SBGGR8),
	PIXFMT(SGBRG8),
	PIXFMT(SGRBG8),
	PIXFMT(SRGGB8),
	PIXFMT(SBGGR10),
	PIXFMT(SGBRG10),
	PIXFMT(SGRBG10),
	PIXFMT(SRGGB10),
	PIXFMT(SBGGR12),
	PIXFMT(SGBRG12),
	PIXFMT(SGRBG12),
	PIXFMT(SRGGB12),
	PIXFMT(SGRBG10DPCM8),
	PIXFMT(SBGGR16),
	PIXFMT(MJPEG),
	PIXFMT(JPEG),
	PIXFMT(DV),
	PIXFMT(MPEG),
	PIXFMT(CPIA1),
	PIXFMT(WNVA),
	PIXFMT(SN9C10X),
	PIXFMT(SN9C20X_I420),
	PIXFMT(PWC1),
	PIXFMT(PWC2),
	PIXFMT(ET61X251),
	PIXFMT(SPCA501),
	PIXFMT(SPCA505),
	PIXFMT(SPCA508),
	PIXFMT(SPCA561),
	PIXFMT(PAC207),
	PIXFMT(MR97310A),
	PIXFMT(SN9C2028),
	PIXFMT(SQ905C),
	PIXFMT(PJPG),
	PIXFMT(OV511),
	PIXFMT(OV518),
	PIXFMT(STV0680),
	PIXFMT(TM6000),
	SYMBOL_END
};

static void print(int lvl, char *msg, ...)
{
	va_list ap;

	if (vars.verbosity < lvl)
		return;

	va_start(ap, msg);
	vprintf(msg, ap);
	va_end(ap);
}

static void error(char *msg, ...)
{
	FILE *f = stdout;
	va_list ap;
	int e = errno;
 
	va_start(ap, msg);
	fprintf(f, "%s: ", name);
	vfprintf(f, msg, ap);
	if (e)
		fprintf(f, ": %s (%i)", strerror(e), e);
	fprintf(f, "\n");
	va_end(ap);
	exit(1);
}

#define xioctl(io, arg) xioctl_(#io, io, arg)

static void xioctl_(char *ios, int ion, void *arg)
{
	int r = ioctl(vars.fd, ion, arg);
	if (r)
		error("%s failed", ios);
}

static void usage(void)
{
	print(1,"Usage: %s [-h] [-d device] [--idsensor] [--idflash]\n", name);
	print(1,"-h	Show this help\n"
		"--help\n"
		"-d	/dev/videoX device node\n"
		"--device\n"
		"-i	Set/get input device (VIDIOC_S/G_INPUT)\n"
		"--input\n"
		"--parm		Set/get parameters (VIDIOC_S/G_PARM)\n"
		"-t		Try format (VIDIOC_TRY_FORMAT)\n"
		"--try_fmt\n"
		"-f		Set/get format (VIDIOC_S/G_FMT)\n"
		"--fmt\n"
		"-r		Request buffers (VIDIOC_REQBUFS)\n"
		"--reqbufs\n"
		"--idsensor	Get sensor identification\n"
		"--idflash	Get flash identification\n"
		"--ctrl-list	List supported V4L2 controls\n"
		"--ctrl=$	Request given V4L2 controls\n"
		"--fmt-list	List supported pixel formats\n"
		"-c $	(VIDIOC_QUERYCAP / VIDIOC_S/G_CTRL / VIDIOC_S/G_EXT_CTRLS)\n"
		"\n"
		"List of V4L2 controls syntax: <[V4L2_CID_]control_name_or_id>[+][=value|?|#][,...]\n"
		"where control_name_or_id is either symbolic name or numerical id.\n"
		"When + is given, use extended controls, otherwise use old-style control call.\n"
		"\"=\" sets the value, \"?\" gets the current value, and \"#\" shows control info.\n");
}

static void symbol_dump(const char *prefix, const struct symbol_list *list)
{
	int i;
	for (i = 0; list[i].symbol != NULL; i++)
		print(0, "%s%s [0x%08X]\n", prefix, list[i].symbol, list[i].id);
}

static int symbol_get(const struct symbol_list *list, const char **symbol)
{
	const char *start;
	const char *end;
	int r, i;

	if (!symbol || !symbol[0])
		error("symbol missing");

	start = *symbol;

	if (isdigit(start[0])) {
		r = strtol(start, (char **)&end, 0);
		if (end == start)
			error("zero-length symbol value");
	} else {
		end = start;
		while (isalnum(*end) || *end == '_') end++;
		if (start == end)
			error("zero-length symbol");
		if (!list)
			error("only numeric value allowed");
		for (i=0; list[i].symbol!=NULL; i++) {
			if (strncasecmp(start, list[i].symbol, end-start) == 0)
				break;
		}
		if (list[i].symbol == NULL)
			error("symbol `%s' not found", start);
		r = list[i].id;
	}

	*symbol = end;
	return r;
}

static void value_get(int *val, int n, const char **ptr)
{
	const char *s = *ptr;
	bool paren = FALSE;
	int i;

	if (!s)
		error("missing integer values");

	if (*s == '(') {
		s++;
		paren = TRUE;
	}

	for (i=0; i<n; i++) {
		char *a;
		val[i] = strtol(s, &a, 0);
		if (a == s)
			error("failure parsing %ith integer argument", i);
		s = a;
		if (i+1 < n) {
			if (*s!=',' && *s!='/' && *s!=';')
				error("missing %ith integer argument", i + 1);
			s++;
		}
	}

	if (paren) {
		if (*s != ')')
			error("missing closing parenthesis");
		s++;
	}

	*ptr = s;
}

static int token_get(const struct token_list *list, const char **token, int val[4])
{
	const char *start;
	const char *end;
	int r, i;

	if (!token || !token[0])
		error("token missing");

	CLEAR(*val);
	end = start = *token;
	while (isalpha(*end) || *end == '.') end++;
	if (start == end)
		error("zero-length token");
	for (i=0; list[i].token!=NULL; i++) {
		if (strncasecmp(start, list[i].token, end-start) == 0)
			break;
	}
	if (list[i].token == NULL)
		error("unrecognized token `%s'", start);
	r = list[i].id;

	if (*end == '=' || *end == ':') {
		if (!(list[i].flags & (TOKEN_F_ARG | TOKEN_F_OPTARG)))
			error("token has argument but should not have");
		end++;
		if (list[i].symbols) {
			val[0] = symbol_get(list[i].symbols, &end);
		} else if (list[i].flags & TOKEN_F_ARG4) {
			value_get(val, 4, &end);
		} else if (list[i].flags & TOKEN_F_ARG2) {
			value_get(val, 2, &end);
		} else {
			value_get(val, 1, &end);
		}
	} else {
		if (list[i].flags & TOKEN_F_ARG)
			error("token does not have an argument but should have");
	}

	/* Go to the beginning of next token in the list */
	if (*end && *end==',') end++;
	*token = end;
	return r;
}

static const char *symbol_str(int id, const struct symbol_list list[])
{
	static char buffer[200];
	int i;

	for (i=0; list[i].symbol; i++)
		if (list[i].id == id)
			break;

	if (list[i].symbol) {
		if (id < 1000)
			sprintf(buffer, "%s [%i]", list[i].symbol, id);
		else
			sprintf(buffer, "%s [0x%08X]", list[i].symbol, id);
	} else {
		sprintf(buffer, "%i", id);
	}
	return buffer;
}

static void vidioc_parm(const char *s)
{
	static const struct symbol_list capturemode[] = {
		{ V4L2_MODE_HIGHQUALITY, "V4L2_MODE_HIGHQUALITY" },
		{ CI_MODE_PREVIEW, "CI_MODE_PREVIEW" },
		{ CI_MODE_VIDEO, "CI_MODE_VIDEO" },
		{ CI_MODE_STILL_CAPTURE, "CI_MODE_STILL_CAPTURE" },
		{ CI_MODE_NONE, "CI_MODE_NONE" },
		SYMBOL_END
	};
	static const struct symbol_list capability[] = {
		{ V4L2_CAP_TIMEPERFRAME, "V4L2_CAP_TIMEPERFRAME" },
		SYMBOL_END
	};
	static const struct token_list list[] = {
		{ 't', TOKEN_F_ARG, "type", buf_types },
		{ 'b', TOKEN_F_ARG, "capability", capability },
		{ 'c', TOKEN_F_ARG, "capturemode", capturemode },
		{ 'f', TOKEN_F_ARG|TOKEN_F_ARG2, "timeperframe", NULL },
		{ 'e', TOKEN_F_ARG, "extendedmode", NULL },
		{ 'r', TOKEN_F_ARG, "readbuffers", NULL },
		{ 'o', TOKEN_F_ARG, "outputmode", capturemode },
		TOKEN_END
	};
	struct v4l2_streamparm p;

	CLEAR(p);
	p.type = vars.reqbufs.type;

	while (*s && *s!='?') {
		int val[4];
		switch (token_get(list, &s, val)) {
		case 't': p.type = val[0]; break;
		case 'b': if (p.type == V4L2_BUF_TYPE_VIDEO_CAPTURE) p.parm.capture.capability = val[0];
			  else  p.parm.output.capability = val[0];
			break;
		case 'c': p.parm.capture.capturemode = val[0]; break;
		case 'f': if (p.type == V4L2_BUF_TYPE_VIDEO_CAPTURE) {
				p.parm.capture.timeperframe.numerator = val[0];
				p.parm.capture.timeperframe.denominator = val[1];
			} else {
				p.parm.output.timeperframe.numerator = val[0];
				p.parm.output.timeperframe.denominator = val[1];
			}
			break;
		case 'e': p.parm.capture.extendedmode = val[0]; break;
		case 'r': p.parm.capture.readbuffers = val[0]; break;
		case 'o': p.parm.output.outputmode = val[0]; break;
		}
	}

	if (*s == '?') {
		print(1, "VIDIOC_G_PARM\n");
		xioctl(VIDIOC_G_PARM, &p);
	} else {
		print(1, "VIDIOC_S_PARM\n");
	}

	print(2, ": type:          %s\n", symbol_str(p.type, buf_types));
	if (p.type == V4L2_BUF_TYPE_VIDEO_CAPTURE) {
		print(2, ": capability:    %s\n", symbol_str(p.parm.capture.capability, capability));
		print(2, ": capturemode:   %s\n", symbol_str(p.parm.capture.capturemode, capturemode));
		print(2, ": timeperframe:  %i/%i\n", p.parm.capture.timeperframe.numerator, p.parm.capture.timeperframe.denominator);
		print(2, ": extendedmode:  %i\n", p.parm.capture.extendedmode);
		print(2, ": readbuffers:   %i\n", p.parm.capture.readbuffers);
	} else if (p.type == V4L2_BUF_TYPE_VIDEO_OUTPUT) {
		print(2, ": capability:    %s\n", symbol_str(p.parm.output.capability, capability));
		print(2, ": outputmode:    %s\n",  symbol_str(p.parm.capture.capturemode, capturemode));
		print(2, ": timeperframe:  %i/%i\n",  p.parm.capture.timeperframe.numerator, p.parm.capture.timeperframe.denominator);
	}

	if (*s != '?') {
		xioctl(VIDIOC_S_PARM, &p);
	}
}

static void vidioc_fmt(bool try, const char *s)
{
	static const struct token_list list[] = {
		{ 't', TOKEN_F_ARG, "type", buf_types },
		{ 'w', TOKEN_F_ARG, "width", NULL },
		{ 'h', TOKEN_F_ARG, "height", NULL },
		{ 'p', TOKEN_F_ARG, "pixelformat", pixelformats },
		{ 'f', TOKEN_F_ARG, "field", NULL },
		{ 'b', TOKEN_F_ARG, "bytesperline", NULL },
		{ 's', TOKEN_F_ARG, "sizeimage", NULL },
		{ 'c', TOKEN_F_ARG, "colorspace", NULL },
		{ 'r', TOKEN_F_ARG, "priv", NULL },
		TOKEN_END
	};
	struct v4l2_format p;

	CLEAR(p);
	p.type = vars.reqbufs.type;

	while (*s && *s!='?') {
		int val[4];
		switch (token_get(list, &s, val)) {
		case 't': p.type = val[0]; break;
		case 'w': p.fmt.pix.width = val[0]; break;
		case 'h': p.fmt.pix.height = val[0]; break;
		case 'p': p.fmt.pix.pixelformat = val[0]; break;
		case 'f': p.fmt.pix.field = val[0]; break;
		case 'b': p.fmt.pix.bytesperline = val[0]; break;
		case 's': p.fmt.pix.sizeimage = val[0]; break;
		case 'c': p.fmt.pix.colorspace = val[0]; break;
		case 'r': p.fmt.pix.priv = val[0]; break;
		}
	}

	if (try) {
		print(1, "VIDIOC_TRY_FMT\n");
		xioctl(VIDIOC_TRY_FMT, &p);
	} else if (*s == '?') {
		print(1, "VIDIOC_G_FMT\n");
		xioctl(VIDIOC_G_FMT, &p);
	} else {
		print(1, "VIDIOC_S_FMT\n");
	}

	print(2, ": type:          %s\n", symbol_str(p.type, buf_types));
	if (p.type == V4L2_BUF_TYPE_VIDEO_CAPTURE) {
		print(2, ": width:         %i\n", p.fmt.pix.width);
		print(2, ": height:        %i\n", p.fmt.pix.height);
		print(2, ": pixelformat:   %s\n", symbol_str(p.fmt.pix.pixelformat, pixelformats));
		print(2, ": field:         %i\n", p.fmt.pix.field);
		print(2, ": bytesperline:  %i\n", p.fmt.pix.bytesperline);
		print(2, ": sizeimage:     %i\n", p.fmt.pix.sizeimage);
		print(2, ": colorspace:    %i\n", p.fmt.pix.colorspace);
		print(2, ": priv:          %i\n", p.fmt.pix.priv);
	}

	if (!try && *s != '?') {
		xioctl(VIDIOC_S_FMT, &p);
	}
}

static void vidioc_reqbufs(const char *s)
{
	static const struct symbol_list memory[] = {
		{ V4L2_MEMORY_MMAP, "MMAP" },
		{ V4L2_MEMORY_USERPTR, "USERPTR" },
		SYMBOL_END
	};
	static const struct token_list list[] = {
		{ 'c', TOKEN_F_ARG, "count", NULL },
		{ 't', TOKEN_F_ARG, "type", buf_types },
		{ 'm', TOKEN_F_ARG, "memory", memory },
		TOKEN_END
	};
	struct v4l2_requestbuffers *p = &vars.reqbufs;

	while (*s && *s!='?') {
		int val[4];
		switch (token_get(list, &s, val)) {
		case 'c': p->count = val[0]; break;
		case 't': p->type = val[0]; break;
		case 'm': p->memory = val[0]; break;
		}
	}

	print(1, "VIDIOC_REQBUFS\n");
	print(2, ": count:   %i\n", p->count);
	print(2, ": type:    %s\n", symbol_str(p->type, buf_types));
	print(2, ": memory:  %s\n", symbol_str(p->memory, memory));
	xioctl(VIDIOC_REQBUFS, p);
}

static __u32 get_control_id(char *name)
{
	__u32 id;

	if (isdigit(*name)) {
		int v;
		if (sscanf(name, "%i", &v) != 1)
			error("bad numeric id");
		id = v;
	} else {
		int i;
		for (i = 0; controls[i].symbol != NULL; i++) {
			if (strcmp(name, controls[i].symbol) == 0)
				break;
			if ((strlen(name) >= sizeof(V4L2_CID)) &&
			    (memcmp(name, V4L2_CID, sizeof(V4L2_CID) - 1) == 0) &&
			    (strcmp(name + sizeof(V4L2_CID) - 1, controls[i].symbol) == 0))
				break;
		}
		if (controls[i].symbol == NULL)
			error("unknown control");
	}

	return id;
}

static const char *get_control_name(__u32 id)
{
	static char buf[11];
	int i;

	for (i = 0; controls[i].symbol != NULL; i++)
		if (controls[i].id == id)
			return controls[i].symbol;

	sprintf(buf, "0x%08X", id);
	return buf;
}

static void close_device()
{
	if (vars.fd == -1)
		return;

	close(vars.fd);
	vars.fd = -1;
	print(1, "CLOSED video device\n");
}

static void open_device(const char *device)
{
	static const char DEFAULT_DEV[] = "/dev/video0";
	if (device == NULL && vars.fd != -1)
		return;

	close_device();	
	if (!device || device[0] == 0)
		device = DEFAULT_DEV;
	print(1, "OPEN video device `%s'\n", device);
	vars.fd = open(device, 0);
	if (vars.fd == -1)
		error("failed to open %s", device);
}

static void v4l2_s_ctrl(__u32 id, __s32 val)
{
	struct v4l2_control c;

	CLEAR(c);
	c.id = id;
	c.value = val;
	print(1, "VIDIOC_S_CTRL[%s] = %i\n", get_control_name(id), c.value);
	xioctl(VIDIOC_S_CTRL, &c);
}

static __s32 v4l2_g_ctrl(__u32 id)
{
	struct v4l2_control c;

	CLEAR(c);
	c.id = id;
	xioctl(VIDIOC_G_CTRL, &c);
	print(1, "VIDIOC_G_CTRL[%s] = %i\n", get_control_name(id), c.value);
	return c.value;
}

static void v4l2_s_ext_ctrl(__u32 id, __s32 val)
{
	struct v4l2_ext_controls cs;
	struct v4l2_ext_control c;

	CLEAR(cs);
	cs.ctrl_class = V4L2_CTRL_ID2CLASS(id);
	cs.count = 1;
	cs.controls = &c;

	CLEAR(c);
	c.id = id;
	c.value = val;

	print(1, "VIDIOC_S_EXT_CTRLS[%s] = %i\n", get_control_name(id), c.value);
	xioctl(VIDIOC_S_EXT_CTRLS, &cs);
}

static void v4l2_query_ctrl(__u32 id)
{
	struct v4l2_queryctrl q;

	CLEAR(q);
	q.id = id;
	xioctl(VIDIOC_QUERYCTRL, &q);
	print(1, "VIDIOC_QUERYCTRL[%s] =\n", get_control_name(id));
	print(2, "  type:    %i\n", q.type);
	print(2, "  name:    %32s\n", q.name);
	print(2, "  limits:  %i..%i / %i\n", q.minimum, q.maximum, q.step);
	print(2, "  default: %i\n", q.default_value);
	print(2, "  flags:   %i\n", q.flags);
}

static __s32 v4l2_g_ext_ctrl(__u32 id)
{
	struct v4l2_ext_controls cs;
	struct v4l2_ext_control c;

	CLEAR(cs);
	cs.ctrl_class = V4L2_CTRL_ID2CLASS(id);
	cs.count = 1;
	cs.controls = &c;

	CLEAR(c);
	c.id = id;

	xioctl(VIDIOC_G_EXT_CTRLS, &cs);
	print(1, "VIDIOC_G_EXT_CTRLS[%s] = %i\n", get_control_name(id), c.value);
	return c.value;
}

static int isident(int c)
{
	return isalnum(c) || c == '_';
}

static void request_controls(char *start)
{
	char *end, *value;
	bool ext, next;
	char op;
	__u32 id;
	int val;

	do {
		for (end = start; isident(*end); end++);
		value = end;
		ext = FALSE;
		if (*value == '+') {
			value++;
			ext = TRUE;
		}
		op = *value++;
		*end = 0;
		next = FALSE;
		id = get_control_id(start);
		if (op == '=') {
			/* Set value */
			for (end = value; isident(*end); end++);
			if (*end == ',')
				next = TRUE;
			if (*end) *end++ = 0;
			if (sscanf(value, "%i", &val) != 1)
				error("bad control value");
			if (ext)
				v4l2_s_ext_ctrl(id, val);
			else
				v4l2_s_ctrl(id, val);
		} else if (op == '?') {
			/* Get value */
			if (*value == ',')
				next = TRUE;
			end = value + 1;
			if (ext)
				v4l2_g_ext_ctrl(id);
			else
				v4l2_g_ctrl(id);
		} else if (op == '#') {
			/* Query control */
			if (*value == ',')
				next = TRUE;
			end = value + 1;
			v4l2_query_ctrl(id);
		} else error("bad request for control");
		start = end;
	} while (next);
}

static void process_options(int argc, char *argv[])
{
	while (1) {
		static const struct option long_options[] = {
			{ "help", 0, NULL, 'h' },
			{ "verbose", 2, NULL, 'v' },
			{ "quiet", 0, NULL, 'q' },
			{ "device", 1, NULL, 'd' },
			{ "input", 1, NULL, 'i' },
			{ "parm", 1, NULL, 'p' },
			{ "try_fmt", 1, NULL, 't' },
			{ "fmt", 1, NULL, 'f' },
			{ "reqbufs", 1, NULL, 'r' },
			{ "idsensor", 0, NULL, 1001 },
			{ "idflash", 0, NULL, 1002 },
			{ "ctrl-list", 0, NULL, 1003 },
			{ "fmt-list", 0, NULL, 1004 },
			{ "ctrl", 1, NULL, 'c' },
			{ 0, 0, 0, 0 }
		};

		int c = getopt_long(argc, argv, "hv::qd:i:p:t:f:r:c:", long_options, NULL);
		if (c == -1)
			break;

		switch (c) {
		case 'h':	/* --help, -h */
			usage();
			return;

		case 'v':	/* --verbose, -v */
			if (optarg) {
				vars.verbosity = atoi(optarg);
			} else {
				vars.verbosity++;
			}
			break;
		
		case 'q':	/* --quiet, -q */
			vars.verbosity--;
			break;

		case 'd':	/* --device, -d */
			open_device(optarg);
			break;

		case 'i': {	/* --input, -i, VIDIOC_S/G_INPUT */
			int i;
			open_device(NULL);
			if (strchr(optarg, '?')) {
				/* G_INPUT */
				xioctl(VIDIOC_G_INPUT, &i);
				print(1, "VIDIOC_G_INPUT -> %i\n", i);
			} else {
				/* S_INPUT */
				i = atoi(optarg);
				print(1, "VIDIOC_S_INPUT <- %i\n", i);
				xioctl(VIDIOC_S_INPUT, &i);
			}
			break;
		}

		case 'p':	/* --parm, -p, VIDIOC_S/G_PARM */
			open_device(NULL);
			vidioc_parm(optarg);
			break;

		case 't':	/* --try-fmt, -t, VIDIOC_TRY_FMT */
			open_device(NULL);
			vidioc_fmt(TRUE, optarg);
			break;

		case 'f':	/* --fmt, -f, VIDIOC_S/G_FMT */
			open_device(NULL);
			vidioc_fmt(FALSE, optarg);
			break;

		case 'r':	/* --reqbufs, -r, VIDIOC_REQBUFS */
			open_device(NULL);
			vidioc_reqbufs(optarg);
			break;

		case 1001: {	/* --idsensor */
			struct atomisp_model_id id;
			open_device(NULL);
			xioctl(ATOMISP_IOC_G_SENSOR_MODEL_ID, &id);
			print(1, "ATOMISP_IOC_G_SENSOR_MODEL_ID: [%s]\n", id.model);
			break;
		}

		case 1002: {	/* --idflash */
			struct atomisp_model_id id;
			open_device(NULL);
			xioctl(ATOMISP_IOC_G_FLASH_MODEL_ID, &id);
			print(1, "ATOMISP_IOC_G_FLASH_MODEL_ID: [%s]\n", id.model);
			break;
		}

		case 1003:	/* --ctrl-list */
			symbol_dump(V4L2_CID, controls);
			return;

		case 1004:	/* --fmt-list */
			symbol_dump(V4L2_PIX_FMT, pixelformats);
			return;

		case 'c':	/* --ctrl, -c, VIDIOC_QUERYCAP / VIDIOC_S/G_CTRL / VIDIOC_S/G_EXT_CTRLS */
			open_device(NULL);
			request_controls(optarg);
			break;

		default:
			error("unknown option");
		}
	}
}

int main(int argc, char *argv[])
{
	print(1, "Starting %s\n", name);
	name = argv[0];
	process_options(argc, argv);
	close_device();
	return 0;
}
