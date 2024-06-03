#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![feature(extern_types)]
extern "C" {
    pub type _IO_wide_data;
    pub type _IO_codecvt;
    pub type _IO_marker;
    fn backtrace(__array: *mut *mut libc::c_void, __size: libc::c_int) -> libc::c_int;
    fn backtrace_symbols_fd(
        __array: *const *mut libc::c_void,
        __size: libc::c_int,
        __fd: libc::c_int,
    );
    static mut optarg: *mut libc::c_char;
    static mut optind: libc::c_int;
    fn getopt_long(
        ___argc: libc::c_int,
        ___argv: *const *mut libc::c_char,
        __shortopts: *const libc::c_char,
        __longopts: *const option,
        __longind: *mut libc::c_int,
    ) -> libc::c_int;
    fn signal(__sig: libc::c_int, __handler: __sighandler_t) -> __sighandler_t;
    fn strtod(_: *const libc::c_char, _: *mut *mut libc::c_char) -> libc::c_double;
    fn strtol(
        _: *const libc::c_char,
        _: *mut *mut libc::c_char,
        _: libc::c_int,
    ) -> libc::c_long;
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
    fn abort() -> !;
    fn exit(_: libc::c_int) -> !;
    fn strrchr(_: *const libc::c_char, _: libc::c_int) -> *mut libc::c_char;
    fn strlen(_: *const libc::c_char) -> libc::c_ulong;
    fn isatty(__fd: libc::c_int) -> libc::c_int;
    static mut stdin: *mut FILE;
    static mut stdout: *mut FILE;
    static mut stderr: *mut FILE;
    fn fclose(__stream: *mut FILE) -> libc::c_int;
    fn fopen(_: *const libc::c_char, _: *const libc::c_char) -> *mut FILE;
    fn fprintf(_: *mut FILE, _: *const libc::c_char, _: ...) -> libc::c_int;
    fn sprintf(_: *mut libc::c_char, _: *const libc::c_char, _: ...) -> libc::c_int;
    fn starcode(
        inputf1: *mut FILE,
        inputf2: *mut FILE,
        outputf1: *mut FILE,
        outputf2: *mut FILE,
        tau: libc::c_int,
        verbose: libc::c_int,
        thrmax: libc::c_int,
        clusteralg: libc::c_int,
        parent_to_child: libc::c_double,
        showclusters: libc::c_int,
        showids: libc::c_int,
        outputt: libc::c_int,
    ) -> libc::c_int;
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct option {
    pub name: *const libc::c_char,
    pub has_arg: libc::c_int,
    pub flag: *mut libc::c_int,
    pub val: libc::c_int,
}
pub type __off_t = libc::c_long;
pub type __off64_t = libc::c_long;
pub type __sighandler_t = Option::<unsafe extern "C" fn(libc::c_int) -> ()>;
pub type size_t = libc::c_ulong;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _IO_FILE {
    pub _flags: libc::c_int,
    pub _IO_read_ptr: *mut libc::c_char,
    pub _IO_read_end: *mut libc::c_char,
    pub _IO_read_base: *mut libc::c_char,
    pub _IO_write_base: *mut libc::c_char,
    pub _IO_write_ptr: *mut libc::c_char,
    pub _IO_write_end: *mut libc::c_char,
    pub _IO_buf_base: *mut libc::c_char,
    pub _IO_buf_end: *mut libc::c_char,
    pub _IO_save_base: *mut libc::c_char,
    pub _IO_backup_base: *mut libc::c_char,
    pub _IO_save_end: *mut libc::c_char,
    pub _markers: *mut _IO_marker,
    pub _chain: *mut _IO_FILE,
    pub _fileno: libc::c_int,
    pub _flags2: libc::c_int,
    pub _old_offset: __off_t,
    pub _cur_column: libc::c_ushort,
    pub _vtable_offset: libc::c_schar,
    pub _shortbuf: [libc::c_char; 1],
    pub _lock: *mut libc::c_void,
    pub _offset: __off64_t,
    pub _codecvt: *mut _IO_codecvt,
    pub _wide_data: *mut _IO_wide_data,
    pub _freeres_list: *mut _IO_FILE,
    pub _freeres_buf: *mut libc::c_void,
    pub __pad5: size_t,
    pub _mode: libc::c_int,
    pub _unused2: [libc::c_char; 20],
}
pub type _IO_lock_t = ();
pub type FILE = _IO_FILE;
pub type C2RustUnnamed = libc::c_uint;
pub const TIDY_OUTPUT: C2RustUnnamed = 3;
pub const NRED_OUTPUT: C2RustUnnamed = 2;
pub const CLUSTER_OUTPUT: C2RustUnnamed = 1;
pub const DEFAULT_OUTPUT: C2RustUnnamed = 0;
pub type C2RustUnnamed_0 = libc::c_uint;
pub const COMPONENTS_CLUSTER: C2RustUnnamed_0 = 2;
pub const SPHERES_CLUSTER: C2RustUnnamed_0 = 1;
pub const MP_CLUSTER: C2RustUnnamed_0 = 0;
#[inline]
unsafe extern "C" fn atoi(mut __nptr: *const libc::c_char) -> libc::c_int {
    return strtol(
        __nptr,
        0 as *mut libc::c_void as *mut *mut libc::c_char,
        10 as libc::c_int,
    ) as libc::c_int;
}
#[inline]
unsafe extern "C" fn atof(mut __nptr: *const libc::c_char) -> libc::c_double {
    return strtod(__nptr, 0 as *mut libc::c_void as *mut *mut libc::c_char);
}
#[no_mangle]
pub static mut USAGE: *mut libc::c_char = b"\nUsage:  starcode [options]\n\n  general options:\n    -d --dist: maximum Levenshtein distance (default auto)\n    -t --threads: number of concurrent threads (default 1)\n    -q --quiet: quiet output (default verbose)\n    -v --version: display version and exit\n\n  cluster options: (default algorithm: message passing)\n    -r --cluster-ratio: min size ratio for merging clusters in\n               message passing (default 5.0)\n    -s --sphere: use sphere clustering algorithm\n    -c --connected-comp: cluster connected components\n\n  input/output options (single file, default)\n    -i --input: input file (default stdin)\n    -o --output: output file (default stdout)\n\n  input options (paired-end fastq files)\n    -1 --input1: input file 1\n    -2 --input2: input file 2\n\n  output options (paired-end fastq files, --non-redundant only)\n       --output1: output file1 (default input1-starcode.fastq)\n       --output2: output file2 (default input2-starcode.fastq)\n\n  output format options\n       --non-redundant: remove redundant sequences from input file(s)\n       --print-clusters: outputs cluster compositions\n       --seq-id: print sequence id numbers (1-based)\n       --tidy: print each sequence and its centroid\n\0"
    as *const u8 as *const libc::c_char as *mut libc::c_char;
#[no_mangle]
pub unsafe extern "C" fn say_usage() {
    fprintf(stderr, b"%s\n\0" as *const u8 as *const libc::c_char, USAGE);
}
#[no_mangle]
pub unsafe extern "C" fn say_version() {
    fprintf(stderr, b"starcode-v1.4\n\0" as *const u8 as *const libc::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn SIGSEGV_handler(mut sig: libc::c_int) {
    let mut array: [*mut libc::c_void; 10] = [0 as *mut libc::c_void; 10];
    let mut size: size_t = 0;
    size = backtrace(array.as_mut_ptr(), 10 as libc::c_int) as size_t;
    fprintf(stderr, b"Error: signal %d:\n\0" as *const u8 as *const libc::c_char, sig);
    backtrace_symbols_fd(array.as_mut_ptr(), size as libc::c_int, 2 as libc::c_int);
    exit(1 as libc::c_int);
}
#[no_mangle]
pub unsafe extern "C" fn outname(mut path: *mut libc::c_char) -> *mut libc::c_char {
    let mut name: *mut libc::c_char = calloc(
        320 as libc::c_int as libc::c_ulong,
        1 as libc::c_int as libc::c_ulong,
    ) as *mut libc::c_char;
    if strlen(path) > 310 as libc::c_int as libc::c_ulong {
        fprintf(
            stderr,
            b"input file name too long (%s)\n\0" as *const u8 as *const libc::c_char,
            path,
        );
        abort();
    }
    let mut c: *mut libc::c_char = strrchr(path, '.' as i32);
    if c.is_null() {
        sprintf(name, b"%s-starcode\0" as *const u8 as *const libc::c_char, path);
    } else {
        *c = '\0' as i32 as libc::c_char;
        sprintf(
            name,
            b"%s-starcode.%s\0" as *const u8 as *const libc::c_char,
            path,
            c.offset(1 as libc::c_int as isize),
        );
        *c = '.' as i32 as libc::c_char;
    }
    return name;
}
unsafe fn main_0(
    mut argc: libc::c_int,
    mut argv: *mut *mut libc::c_char,
) -> libc::c_int {
    signal(
        11 as libc::c_int,
        Some(SIGSEGV_handler as unsafe extern "C" fn(libc::c_int) -> ()),
    );
    static mut nr_flag: libc::c_int = 0 as libc::c_int;
    static mut td_flag: libc::c_int = 0 as libc::c_int;
    static mut sp_flag: libc::c_int = 0 as libc::c_int;
    static mut vb_flag: libc::c_int = 1 as libc::c_int;
    static mut cl_flag: libc::c_int = 0 as libc::c_int;
    static mut id_flag: libc::c_int = 0 as libc::c_int;
    static mut cp_flag: libc::c_int = 0 as libc::c_int;
    let mut dist: libc::c_int = -(1 as libc::c_int);
    let mut threads: libc::c_int = -(1 as libc::c_int);
    let mut cluster_ratio: libc::c_double = -(1 as libc::c_int) as libc::c_double;
    let UNSET: *mut libc::c_char = b"unset\0" as *const u8 as *const libc::c_char
        as *mut libc::c_char;
    let mut input: *mut libc::c_char = UNSET;
    let mut input1: *mut libc::c_char = UNSET;
    let mut input2: *mut libc::c_char = UNSET;
    let mut output: *mut libc::c_char = UNSET;
    let mut output1: *mut libc::c_char = UNSET;
    let mut output2: *mut libc::c_char = UNSET;
    if argc == 1 as libc::c_int && isatty(0 as libc::c_int) != 0 {
        say_usage();
        return 0 as libc::c_int;
    }
    let mut c: libc::c_int = 0;
    loop {
        let mut option_index: libc::c_int = 0 as libc::c_int;
        static mut long_options: [option; 19] = unsafe {
            [
                {
                    let mut init = option {
                        name: b"print-clusters\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &cl_flag as *const libc::c_int as *mut libc::c_int,
                        val: 1 as libc::c_int,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"seq-id\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &id_flag as *const libc::c_int as *mut libc::c_int,
                        val: 1 as libc::c_int,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"non-redundant\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &nr_flag as *const libc::c_int as *mut libc::c_int,
                        val: 1 as libc::c_int,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"tidy\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &td_flag as *const libc::c_int as *mut libc::c_int,
                        val: 1 as libc::c_int,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"quiet\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &vb_flag as *const libc::c_int as *mut libc::c_int,
                        val: 0 as libc::c_int,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"sphere\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &sp_flag as *const libc::c_int as *mut libc::c_int,
                        val: 's' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"connected-comp\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: &cp_flag as *const libc::c_int as *mut libc::c_int,
                        val: 'c' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"version\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 'v' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"dist\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 'd' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"cluster-ratio\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 'r' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"help\0" as *const u8 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 'h' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"input\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 'i' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"input1\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: '1' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"input2\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: '2' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"output\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 'o' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"threads\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 't' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"output1\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: '3' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: b"output2\0" as *const u8 as *const libc::c_char,
                        has_arg: 1 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: '4' as i32,
                    };
                    init
                },
                {
                    let mut init = option {
                        name: 0 as *const libc::c_char,
                        has_arg: 0 as libc::c_int,
                        flag: 0 as *const libc::c_int as *mut libc::c_int,
                        val: 0 as libc::c_int,
                    };
                    init
                },
            ]
        };
        c = getopt_long(
            argc,
            argv,
            b"1:2:3:4:d:hi:o:qcst:r:v\0" as *const u8 as *const libc::c_char,
            long_options.as_mut_ptr(),
            &mut option_index,
        );
        if c == -(1 as libc::c_int) {
            break;
        }
        match c {
            0 => {}
            49 => {
                if input1 == UNSET {
                    input1 = optarg;
                } else {
                    fprintf(
                        stderr,
                        b"%s --input1 set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            50 => {
                if input2 == UNSET {
                    input2 = optarg;
                } else {
                    fprintf(
                        stderr,
                        b"%s --input2 set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            51 => {
                if output1 == UNSET {
                    output1 = optarg;
                } else {
                    fprintf(
                        stderr,
                        b"%s --output1 set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            52 => {
                if output2 == UNSET {
                    output2 = optarg;
                } else {
                    fprintf(
                        stderr,
                        b"%s --output2 set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            100 => {
                if dist < 0 as libc::c_int {
                    dist = atoi(optarg);
                    if dist > 8 as libc::c_int {
                        fprintf(
                            stderr,
                            b"%s --dist cannot exceed %d\n\0" as *const u8
                                as *const libc::c_char,
                            b"starcode error:\0" as *const u8 as *const libc::c_char,
                            8 as libc::c_int,
                        );
                        return 1 as libc::c_int;
                    }
                } else {
                    fprintf(
                        stderr,
                        b"%s --distance set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            104 => {
                say_version();
                say_usage();
                return 0 as libc::c_int;
            }
            105 => {
                if input == UNSET {
                    input = optarg;
                } else {
                    fprintf(
                        stderr,
                        b"%s --input set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            111 => {
                if output == UNSET {
                    output = optarg;
                } else {
                    fprintf(
                        stderr,
                        b"%s --output set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            113 => {
                vb_flag = 0 as libc::c_int;
            }
            115 => {
                sp_flag = 1 as libc::c_int;
            }
            99 => {
                cp_flag = 1 as libc::c_int;
            }
            116 => {
                if threads < 0 as libc::c_int {
                    threads = atoi(optarg);
                    if threads < 1 as libc::c_int {
                        fprintf(
                            stderr,
                            b"%s --threads must be a positive integer\n\0" as *const u8
                                as *const libc::c_char,
                            b"starcode error:\0" as *const u8 as *const libc::c_char,
                        );
                        say_usage();
                        return 1 as libc::c_int;
                    }
                } else {
                    fprintf(
                        stderr,
                        b"%s --thread set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            114 => {
                if cluster_ratio < 0 as libc::c_int as libc::c_double {
                    cluster_ratio = atof(optarg);
                    if cluster_ratio < 1 as libc::c_int as libc::c_double {
                        fprintf(
                            stderr,
                            b"%s --cluster-ratio must be greater or equal than 1.0.\n\0"
                                as *const u8 as *const libc::c_char,
                            b"starcode error:\0" as *const u8 as *const libc::c_char,
                        );
                        say_usage();
                        return 1 as libc::c_int;
                    }
                } else {
                    fprintf(
                        stderr,
                        b"%s --cluster-ratio set more than once\n\0" as *const u8
                            as *const libc::c_char,
                        b"starcode error:\0" as *const u8 as *const libc::c_char,
                    );
                    say_usage();
                    return 1 as libc::c_int;
                }
            }
            118 => {
                say_version();
                return 0 as libc::c_int;
            }
            _ => {
                say_usage();
                return 1 as libc::c_int;
            }
        }
    }
    if optind < argc {
        if optind == argc - 1 as libc::c_int && (input == UNSET && input1 == UNSET) {
            input = *argv.offset(optind as isize);
        } else {
            fprintf(
                stderr,
                b"%s too many options\n\0" as *const u8 as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
            );
            say_usage();
            return 1 as libc::c_int;
        }
    }
    if nr_flag != 0 && (cl_flag != 0 || id_flag != 0) {
        fprintf(
            stderr,
            b"%s --non-redundant flag is incompatible with --print-clusters and --seq-id\n\0"
                as *const u8 as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    if input != UNSET && (input1 != UNSET || input2 != UNSET) {
        fprintf(
            stderr,
            b"%s --input and --input1/2 are incompatible\n\0" as *const u8
                as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    if input1 == UNSET && input2 != UNSET {
        fprintf(
            stderr,
            b"%s --input2 set without --input1\n\0" as *const u8 as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    if input2 == UNSET && input1 != UNSET {
        fprintf(
            stderr,
            b"%s --input1 set without --input2\n\0" as *const u8 as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    if nr_flag != 0 && output != UNSET && (input1 != UNSET || input2 != UNSET) {
        fprintf(
            stderr,
            b"%s cannot specify --output for paired-end fastq file with --non-redundant\n\0"
                as *const u8 as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    if sp_flag != 0 && cp_flag != 0 {
        fprintf(
            stderr,
            b"%s --sphere and --connected-comp are incompatible\n\0" as *const u8
                as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    if td_flag != 0 && (nr_flag != 0 || cl_flag != 0 || id_flag != 0) {
        fprintf(
            stderr,
            b"%s --tidy flag is not compatible with options --print-clusters, --seq-id and --non-redundant\n\0"
                as *const u8 as *const libc::c_char,
            b"starcode error:\0" as *const u8 as *const libc::c_char,
        );
        say_usage();
        return 1 as libc::c_int;
    }
    let mut output_type: libc::c_int = 0;
    if nr_flag != 0 {
        output_type = NRED_OUTPUT as libc::c_int;
    } else if td_flag != 0 {
        output_type = TIDY_OUTPUT as libc::c_int;
    } else {
        output_type = DEFAULT_OUTPUT as libc::c_int;
    }
    let mut cluster_alg: libc::c_int = 0;
    if cp_flag != 0 {
        cluster_alg = COMPONENTS_CLUSTER as libc::c_int;
    } else if sp_flag != 0 {
        cluster_alg = SPHERES_CLUSTER as libc::c_int;
    } else {
        cluster_alg = MP_CLUSTER as libc::c_int;
    }
    let mut inputf1: *mut FILE = 0 as *mut FILE;
    let mut inputf2: *mut FILE = 0 as *mut FILE;
    let mut outputf1: *mut FILE = 0 as *mut FILE;
    let mut outputf2: *mut FILE = 0 as *mut FILE;
    if input != UNSET {
        inputf1 = fopen(input, b"r\0" as *const u8 as *const libc::c_char);
        if inputf1.is_null() {
            fprintf(
                stderr,
                b"%s cannot open file %s\n\0" as *const u8 as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
                input,
            );
            say_usage();
            return 1 as libc::c_int;
        }
    } else if input1 != UNSET {
        inputf1 = fopen(input1, b"r\0" as *const u8 as *const libc::c_char);
        if inputf1.is_null() {
            fprintf(
                stderr,
                b"%s cannot open file %s\n\0" as *const u8 as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
                input1,
            );
            say_usage();
            return 1 as libc::c_int;
        }
        inputf2 = fopen(input2, b"r\0" as *const u8 as *const libc::c_char);
        if inputf2.is_null() {
            fprintf(
                stderr,
                b"%s cannot open file %s\n\0" as *const u8 as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
                input2,
            );
            say_usage();
            return 1 as libc::c_int;
        }
    } else {
        inputf1 = stdin;
    }
    if output != UNSET {
        outputf1 = fopen(output, b"w\0" as *const u8 as *const libc::c_char);
        if outputf1.is_null() {
            fprintf(
                stderr,
                b"%s cannot write to file %s\n\0" as *const u8 as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
                output,
            );
            say_usage();
            return 1 as libc::c_int;
        }
    } else if nr_flag != 0 && input1 != UNSET && input2 != UNSET {
        if output1 == UNSET {
            output1 = outname(input1);
            outputf1 = fopen(output1, b"w\0" as *const u8 as *const libc::c_char);
            free(output1 as *mut libc::c_void);
            output1 = 0 as *mut libc::c_char;
        } else {
            outputf1 = fopen(output1, b"w\0" as *const u8 as *const libc::c_char);
        }
        if outputf1.is_null() {
            fprintf(
                stderr,
                b"%s cannot write to output file 1\n\0" as *const u8
                    as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
            );
            say_usage();
            return 1 as libc::c_int;
        }
        if output2 == UNSET {
            output2 = outname(input2);
            outputf2 = fopen(output2, b"w\0" as *const u8 as *const libc::c_char);
            free(output2 as *mut libc::c_void);
            output2 = 0 as *mut libc::c_char;
        } else {
            outputf2 = fopen(output2, b"w\0" as *const u8 as *const libc::c_char);
        }
        if outputf2.is_null() {
            fprintf(
                stderr,
                b"%s cannot write to output file 2\n\0" as *const u8
                    as *const libc::c_char,
                b"starcode error:\0" as *const u8 as *const libc::c_char,
            );
            say_usage();
            return 1 as libc::c_int;
        }
    } else {
        outputf1 = stdout;
    }
    if threads < 0 as libc::c_int {
        threads = 1 as libc::c_int;
    }
    if cluster_ratio < 0 as libc::c_int as libc::c_double {
        cluster_ratio = 5 as libc::c_int as libc::c_double;
    }
    if cluster_ratio == 1.0f64 && vb_flag != 0 {
        fprintf(
            stderr,
            b"warning: setting cluster-ratio to 1.0 may result in arbitrary cluster breaks.\n\0"
                as *const u8 as *const libc::c_char,
        );
    }
    let mut exitcode: libc::c_int = starcode(
        inputf1,
        inputf2,
        outputf1,
        outputf2,
        dist,
        vb_flag,
        threads,
        cluster_alg,
        cluster_ratio,
        cl_flag,
        id_flag,
        output_type,
    );
    if inputf1 != stdin {
        fclose(inputf1);
    }
    if !inputf2.is_null() {
        fclose(inputf2);
    }
    if outputf1 != stdout {
        fclose(outputf1);
    }
    if !outputf2.is_null() {
        fclose(outputf2);
    }
    return exitcode;
}
pub fn main() {
    let mut args: Vec::<*mut libc::c_char> = Vec::new();
    for arg in ::std::env::args() {
        args.push(
            (::std::ffi::CString::new(arg))
                .expect("Failed to convert argument into CString.")
                .into_raw(),
        );
    }
    args.push(::core::ptr::null_mut());
    unsafe {
        ::std::process::exit(
            main_0(
                (args.len() - 1) as libc::c_int,
                args.as_mut_ptr() as *mut *mut libc::c_char,
            ) as i32,
        )
    }
}
