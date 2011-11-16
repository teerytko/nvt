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
#include <linux/videodev2.h>

#include "linux/atomisp.h"

char *name = "v4l2n";

#define FALSE		0
#define TRUE		(!FALSE)

#define STRINGIFY_(x)	#x
#define STRINGIFY(x)	STRINGIFY_(x)

#define CLEAR(x)	memset(&(x), 0, sizeof(x));
#define SIZE(x)		(sizeof(x)/sizeof((x)[0]))

typedef unsigned char bool;

static struct {
	int verbosity;
	int fd;
} vars = {
	.verbosity = 2,
	.fd = -1,
};

#define V4L2_CID	"V4L2_CID_"
#define CONTROL(id)	{ V4L2_CID_##id, (#id) }
static struct {
	__u32 id;
	char *name;
} controls[] = {
	CONTROL(BRIGHTNESS),

	/* Flash controls */
	CONTROL(FLASH_DURATION),
	CONTROL(FLASH_INTENSITY),
	CONTROL(TORCH_INTENSITY),
	CONTROL(INDICATOR_INTENSITY),
	CONTROL(FLASH_TRIGGER),
	CONTROL(FLASH_MODE),
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
		"-i	Set/get input device\n"
		"--input\n"
		"--idsensor	Get sensor identification\n"
		"--idflash	Get flash identification\n"
		"--ctrl-list	List supported V4L2 controls\n"
		"--ctrl=$	Request given V4L2 controls\n"
		"-c $\n"
		"\n"
		"List of V4L2 controls syntax: <[V4L2_CID_]control_name_or_id>[+][=value|?|#][,...]\n"
		"where control_name_or_id is either symbolic name or numerical id.\n"
		"When + is given, use extended controls, otherwise use old-style control call.\n"
		"\"=\" sets the value, \"?\" gets the current value, and \"#\" shows control info.\n");
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
		for (i = 0; i < SIZE(controls); i++) {
			if (strcmp(name, controls[i].name) == 0)
				break;
			if ((strlen(name) >= sizeof(V4L2_CID)) &&
			    (memcmp(name, V4L2_CID, sizeof(V4L2_CID) - 1) == 0) &&
			    (strcmp(name + sizeof(V4L2_CID) - 1, controls[i].name) == 0))
				break;
		}
		if (i >= SIZE(controls))
			error("unknown control");
	}

	return id;
}

static char *get_control_name(__u32 id)
{
	static char buf[11];
	int i;

	for (i = 0; i < SIZE(controls); i++)
		if (controls[i].id == id)
			return controls[i].name;

	sprintf(buf, "0x%08X", id);
	return buf;
}

static void close_device()
{
	if (vars.fd != -1)
		close(vars.fd);
	vars.fd = -1;
}

static void open_device(const char *device)
{
	static const char DEFAULT_DEV[] = "/dev/video0";
	if (device == NULL && vars.fd != -1)
		return;

	close_device();	
	if (!device || device[0] == 0)
		device = DEFAULT_DEV;
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
	xioctl(VIDIOC_S_CTRL, &c);
	print(1, "VIDIOC_S_CTRL[%s] = %i\n", get_control_name(id), c.value);
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

	xioctl(VIDIOC_S_EXT_CTRLS, &cs);
	print(1, "VIDIOC_S_EXT_CTRLS[%s] = %i\n", get_control_name(id), c.value);
}

static void v4l2_query_ctrl(__u32 id)
{
	struct v4l2_queryctrl q;

	CLEAR(q);
	q.id = id;
	xioctl(VIDIOC_QUERYCTRL, &q);
	print(1, "VIDIOC_QUERYCTRL[%s] =\n", get_control_name(id));
	print(1, "  type:    %i\n", q.type);
	print(1, "  name:    %32s\n", q.name);
	print(1, "  limits:  %i..%i / %i\n", q.minimum, q.maximum, q.step);
	print(1, "  default: %i\n", q.default_value);
	print(1, "  flags:   %i\n", q.flags);
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
			{ "idsensor", 0, NULL, 1001 },
			{ "idflash", 0, NULL, 1002 },
			{ "ctrl-list", 0, NULL, 1003 },
			{ "ctrl", 1, NULL, 'c' },
			{ 0, 0, 0, 0 }
		};

		int c = getopt_long(argc, argv, "hv::qd:i:c:", long_options, NULL);
		if (c == -1)
			break;

		switch (c) {
		case 'h':
			usage();
			return;

		case 'v':
			if (optarg) {
				vars.verbosity = atoi(optarg);
			} else {
				vars.verbosity++;
			}
			break;
		
		case 'q':
			vars.verbosity--;
			break;

		case 'd':
			open_device(optarg);
			break;

		case 'i': {	/* S/G_INPUT */
			int i;
			open_device(NULL);
			if (optarg[0] == '?') {
				/* G_INPUT */
				xioctl(VIDIOC_G_INPUT, &i);
				print(1, "VIDIOC_G_INPUT -> %i\n", i);
			} else {
				i = atoi(optarg);
				print(1, "VIDIOC_S_INPUT <- %i\n", i);
				xioctl(VIDIOC_S_INPUT, &i);
			}
			break;
		}

		case 1001: {
			struct atomisp_model_id id;
			open_device(NULL);
			xioctl(ATOMISP_IOC_G_SENSOR_MODEL_ID, &id);
			print(1, "ATOMISP_IOC_G_SENSOR_MODEL_ID: [%s]\n", id.model);
			break;
		}

		case 1002: {
			struct atomisp_model_id id;
			open_device(NULL);
			xioctl(ATOMISP_IOC_G_FLASH_MODEL_ID, &id);
			print(1, "ATOMISP_IOC_G_FLASH_MODEL_ID: [%s]\n", id.model);
			break;
		}

		case 1003: {
			int i;
			for (i = 0; i < SIZE(controls); i++)
				print(1, "V4L2_CID_%s [0x%08X]\n", controls[i].name, controls[i].id);
			exit(0);
		}

		case 'c':
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
