#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
extern "C" {



    pub type _IO_wide_data;


    pub type _IO_codecvt;

    pub type _IO_marker;
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn realloc(_: *mut libc::c_void, _: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
    fn memcpy(
        _: *mut libc::c_void,
        _: *const libc::c_void,
        _: libc::c_ulong,
    ) -> *mut libc::c_void;
    fn memset(
        _: *mut libc::c_void,
        _: libc::c_int,
        _: libc::c_ulong,
    ) -> *mut libc::c_void;
    static mut stderr: *mut FILE;
    fn fprintf(_: *mut FILE, _: *const libc::c_char, _: ...) -> libc::c_int;
    fn strlen(_: *const libc::c_char) -> libc::c_ulong;
}
pub type size_t = libc::c_ulong;
pub type __int32_t = libc::c_int;
pub type __uint32_t = libc::c_uint;
pub type __off_t = libc::c_long;
pub type __off64_t = libc::c_long;
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
pub type int32_t = __int32_t;
pub type uint32_t = __uint32_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct gstack_t {
    pub nslots: size_t,
    pub nitems: size_t,
    pub items: [*mut libc::c_void; 0],
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct info_t {
    pub height: libc::c_uint,
    pub pebbles: *mut *mut gstack_t,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct node_t {
    pub child: [*mut libc::c_void; 6],
    pub path: uint32_t,
    pub cache: [libc::c_char; 17],
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct trie_t {
    pub root: *mut node_t,
    pub info: *mut info_t,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct arg_t {
    pub hits: *mut *mut gstack_t,
    pub pebbles: *mut *mut gstack_t,
    pub tau: libc::c_char,
    pub maxtau: libc::c_char,
    pub query: *mut libc::c_int,
    pub seed_depth: libc::c_int,
    pub height: libc::c_int,
    pub err: libc::c_int,
}
static mut translate: [libc::c_int; 256] = [
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    5 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    1 as libc::c_int,
    0 as libc::c_int,
    2 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    3 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    4 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    1 as libc::c_int,
    0 as libc::c_int,
    2 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    3 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    4 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
    0 as libc::c_int,
];
static mut altranslate: [libc::c_int; 256] = [
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    5 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    0 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    1 as libc::c_int,
    6 as libc::c_int,
    2 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    3 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    4 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    1 as libc::c_int,
    6 as libc::c_int,
    2 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    3 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    4 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
    6 as libc::c_int,
];
#[no_mangle]
pub static mut ERROR: libc::c_int = 0 as libc::c_int;
#[no_mangle]
pub static mut TOWER_TOP: *mut gstack_t = 0 as *const gstack_t as *mut gstack_t;
#[no_mangle]
pub unsafe extern "C" fn get_height(mut trie: *mut trie_t) -> libc::c_int {
    return (*(*trie).info).height as libc::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn search(
    mut trie: *mut trie_t,
    mut query: *const libc::c_char,
    tau: libc::c_int,
    mut hits: *mut *mut gstack_t,
    mut start_depth: libc::c_int,
    seed_depth: libc::c_int,
) -> libc::c_int {
    ERROR = 0 as libc::c_int;
    let mut height: libc::c_int = get_height(trie);
    if tau > 8 as libc::c_int {
        fprintf(
            stderr,
            b"error: requested tau greater than %d\n\0" as *const u8
                as *const libc::c_char,
            8 as libc::c_int,
        );
        return 155 as libc::c_int;
    }
    let mut length: libc::c_int = strlen(query) as libc::c_int;
    if length > height {
        fprintf(
            stderr,
            b"error: query longer than allowed max\n\0" as *const u8
                as *const libc::c_char,
        );
        return 161 as libc::c_int;
    }
    let mut info: *mut info_t = (*trie).info;
    start_depth = if start_depth > 0 as libc::c_int {
        start_depth
    } else {
        0 as libc::c_int
    };
    let mut i: libc::c_int = start_depth + 1 as libc::c_int;
    while i <= (if seed_depth < height { seed_depth } else { height }) {
        (**((*info).pebbles).offset(i as isize)).nitems = 0 as libc::c_int as size_t;
        i += 1;
        i;
    }
    let mut translated: [libc::c_int; 1024] = [0; 1024];
    translated[0 as libc::c_int as usize] = length;
    translated[(length + 1 as libc::c_int) as usize] = -(1 as libc::c_int);
    let mut i_0: libc::c_int = if 0 as libc::c_int > start_depth - 8 as libc::c_int {
        0 as libc::c_int
    } else {
        start_depth - 8 as libc::c_int
    };
    while i_0 < length {
        translated[(i_0 + 1 as libc::c_int)
            as usize] = altranslate[*query.offset(i_0 as isize) as libc::c_int as usize];
        i_0 += 1;
        i_0;
    }
    let mut arg: arg_t = {
        let mut init = arg_t {
            hits: hits,
            pebbles: (*info).pebbles,
            tau: tau as libc::c_char,
            maxtau: 0,
            query: translated.as_mut_ptr(),
            seed_depth: seed_depth,
            height: height,
            err: 0,
        };
        init
    };
    let mut pebbles: *mut gstack_t = *((*info).pebbles).offset(start_depth as isize);
    let mut i_1: libc::c_uint = 0 as libc::c_int as libc::c_uint;
    while (i_1 as libc::c_ulong) < (*pebbles).nitems {
        let mut start_node: *mut node_t = *((*pebbles).items)
            .as_mut_ptr()
            .offset(i_1 as isize) as *mut node_t;
        poucet(start_node, start_depth + 1 as libc::c_int, arg);
        i_1 = i_1.wrapping_add(1);
        i_1;
    }
    return check_trie_error_and_reset();
}
#[no_mangle]
pub unsafe extern "C" fn poucet(
    mut node: *mut node_t,
    depth: libc::c_int,
    mut arg: arg_t,
) {
    let mut pcache: *mut libc::c_char = ((*node).cache)
        .as_mut_ptr()
        .offset(8 as libc::c_int as isize);
    let mut maxa: libc::c_int = if (depth - 1 as libc::c_int) < arg.tau as libc::c_int {
        depth - 1 as libc::c_int
    } else {
        arg.tau as libc::c_int
    };
    let mut mmatch: libc::c_uchar = 0;
    let mut shift: libc::c_uchar = 0;
    let mut common: [libc::c_char; 9] = [
        1 as libc::c_int as libc::c_char,
        2 as libc::c_int as libc::c_char,
        3 as libc::c_int as libc::c_char,
        4 as libc::c_int as libc::c_char,
        5 as libc::c_int as libc::c_char,
        6 as libc::c_int as libc::c_char,
        7 as libc::c_int as libc::c_char,
        8 as libc::c_int as libc::c_char,
        9 as libc::c_int as libc::c_char,
    ];
    let mut path: int32_t = (*node).path as int32_t;
    if maxa > 0 as libc::c_int {
        mmatch = ((if *(arg.query).offset((depth - 1 as libc::c_int) as isize)
            == 5 as libc::c_int
        {
            0 as libc::c_int
        } else {
            *pcache.offset(maxa as isize) as libc::c_int
        })
            + (path >> 4 as libc::c_int * (maxa - 1 as libc::c_int) & 15 as libc::c_int
                != *(arg.query).offset(depth as isize)) as libc::c_int) as libc::c_uchar;
        shift = ((if (*pcache.offset((maxa - 1 as libc::c_int) as isize) as libc::c_int)
            < common[maxa as usize] as libc::c_int
        {
            *pcache.offset((maxa - 1 as libc::c_int) as isize) as libc::c_int
        } else {
            common[maxa as usize] as libc::c_int
        }) + 1 as libc::c_int) as libc::c_uchar;
        common[(maxa - 1 as libc::c_int)
            as usize] = (if (mmatch as libc::c_int) < shift as libc::c_int {
            mmatch as libc::c_int
        } else {
            shift as libc::c_int
        }) as libc::c_char;
        let mut a: libc::c_int = maxa - 1 as libc::c_int;
        while a > 0 as libc::c_int {
            mmatch = (*pcache.offset(a as isize) as libc::c_int
                + (path >> 4 as libc::c_int * (a - 1 as libc::c_int) & 15 as libc::c_int
                    != *(arg.query).offset(depth as isize)) as libc::c_int)
                as libc::c_uchar;
            shift = ((if (*pcache.offset((a - 1 as libc::c_int) as isize) as libc::c_int)
                < common[a as usize] as libc::c_int
            {
                *pcache.offset((a - 1 as libc::c_int) as isize) as libc::c_int
            } else {
                common[a as usize] as libc::c_int
            }) + 1 as libc::c_int) as libc::c_uchar;
            common[(a - 1 as libc::c_int)
                as usize] = (if (mmatch as libc::c_int) < shift as libc::c_int {
                mmatch as libc::c_int
            } else {
                shift as libc::c_int
            }) as libc::c_char;
            a -= 1;
            a;
        }
    }
    let mut child: *mut node_t = 0 as *mut node_t;
    let mut current_block_33: u64;
    let mut i: libc::c_int = 0 as libc::c_int;
    while i < 6 as libc::c_int {
        child = (*node).child[i as usize] as *mut node_t;
        if !child.is_null() {
            let mut local_cache: [libc::c_char; 19] = [
                9 as libc::c_int as libc::c_char,
                8 as libc::c_int as libc::c_char,
                7 as libc::c_int as libc::c_char,
                6 as libc::c_int as libc::c_char,
                5 as libc::c_int as libc::c_char,
                4 as libc::c_int as libc::c_char,
                3 as libc::c_int as libc::c_char,
                2 as libc::c_int as libc::c_char,
                1 as libc::c_int as libc::c_char,
                0 as libc::c_int as libc::c_char,
                1 as libc::c_int as libc::c_char,
                2 as libc::c_int as libc::c_char,
                3 as libc::c_int as libc::c_char,
                4 as libc::c_int as libc::c_char,
                5 as libc::c_int as libc::c_char,
                6 as libc::c_int as libc::c_char,
                7 as libc::c_int as libc::c_char,
                8 as libc::c_int as libc::c_char,
                9 as libc::c_int as libc::c_char,
            ];
            let mut ccache: *mut libc::c_char = if depth == arg.height {
                local_cache.as_mut_ptr().offset(9 as libc::c_int as isize)
            } else {
                ((*child).cache).as_mut_ptr().offset(8 as libc::c_int as isize)
            };
            memcpy(
                ccache.offset(1 as libc::c_int as isize) as *mut libc::c_void,
                common.as_mut_ptr() as *const libc::c_void,
                (8 as libc::c_int as libc::c_ulong)
                    .wrapping_mul(
                        ::core::mem::size_of::<libc::c_char>() as libc::c_ulong,
                    ),
            );
            if maxa > 0 as libc::c_int {
                mmatch = ((if path & 15 as libc::c_int == 5 as libc::c_int {
                    0 as libc::c_int
                } else {
                    *pcache.offset(-maxa as isize) as libc::c_int
                }) + (i != *(arg.query).offset((depth - maxa) as isize)) as libc::c_int)
                    as libc::c_uchar;
                shift = ((if (*pcache.offset((1 as libc::c_int - maxa) as isize)
                    as libc::c_int) < maxa + 1 as libc::c_int
                {
                    *pcache.offset((1 as libc::c_int - maxa) as isize) as libc::c_int
                } else {
                    maxa + 1 as libc::c_int
                }) + 1 as libc::c_int) as libc::c_uchar;
                *ccache
                    .offset(
                        -maxa as isize,
                    ) = (if (mmatch as libc::c_int) < shift as libc::c_int {
                    mmatch as libc::c_int
                } else {
                    shift as libc::c_int
                }) as libc::c_char;
                let mut a_0: libc::c_int = maxa - 1 as libc::c_int;
                while a_0 > 0 as libc::c_int {
                    mmatch = (*pcache.offset(-a_0 as isize) as libc::c_int
                        + (i != *(arg.query).offset((depth - a_0) as isize))
                            as libc::c_int) as libc::c_uchar;
                    shift = ((if (*pcache.offset((1 as libc::c_int - a_0) as isize)
                        as libc::c_int)
                        < *ccache.offset((-a_0 - 1 as libc::c_int) as isize)
                            as libc::c_int
                    {
                        *pcache.offset((1 as libc::c_int - a_0) as isize) as libc::c_int
                    } else {
                        *ccache.offset((-a_0 - 1 as libc::c_int) as isize) as libc::c_int
                    }) + 1 as libc::c_int) as libc::c_uchar;
                    *ccache
                        .offset(
                            -a_0 as isize,
                        ) = (if (mmatch as libc::c_int) < shift as libc::c_int {
                        mmatch as libc::c_int
                    } else {
                        shift as libc::c_int
                    }) as libc::c_char;
                    a_0 -= 1;
                    a_0;
                }
            }
            mmatch = (*pcache.offset(0 as libc::c_int as isize) as libc::c_int
                + (i != *(arg.query).offset(depth as isize)) as libc::c_int)
                as libc::c_uchar;
            shift = ((if (*ccache.offset(-(1 as libc::c_int) as isize) as libc::c_int)
                < *ccache.offset(1 as libc::c_int as isize) as libc::c_int
            {
                *ccache.offset(-(1 as libc::c_int) as isize) as libc::c_int
            } else {
                *ccache.offset(1 as libc::c_int as isize) as libc::c_int
            }) + 1 as libc::c_int) as libc::c_uchar;
            *ccache
                .offset(
                    0 as libc::c_int as isize,
                ) = (if (mmatch as libc::c_int) < shift as libc::c_int {
                mmatch as libc::c_int
            } else {
                shift as libc::c_int
            }) as libc::c_char;
            if !(*ccache.offset(0 as libc::c_int as isize) as libc::c_int
                > arg.tau as libc::c_int)
            {
                if depth == arg.height {
                    if push(
                        child as *mut libc::c_void,
                        (arg.hits)
                            .offset(
                                *ccache.offset(0 as libc::c_int as isize) as libc::c_int
                                    as isize,
                            ),
                    ) != 0
                    {
                        ERROR = 335 as libc::c_int;
                    }
                } else {
                    if depth <= arg.seed_depth {
                        if push(
                            child as *mut libc::c_void,
                            (arg.pebbles).offset(depth as isize),
                        ) != 0
                        {
                            ERROR = 341 as libc::c_int;
                        }
                    }
                    if depth > arg.seed_depth {
                        let mut can_dash: libc::c_int = 1 as libc::c_int;
                        let mut a_1: libc::c_int = -maxa;
                        while a_1 < maxa + 1 as libc::c_int {
                            if (*ccache.offset(a_1 as isize) as libc::c_int)
                                < arg.tau as libc::c_int
                            {
                                can_dash = 0 as libc::c_int;
                                break;
                            } else {
                                a_1 += 1;
                                a_1;
                            }
                        }
                        if can_dash != 0 {
                            dash(
                                child,
                                (arg.query)
                                    .offset(depth as isize)
                                    .offset(1 as libc::c_int as isize),
                                arg,
                            );
                            current_block_33 = 2968425633554183086;
                        } else {
                            current_block_33 = 3934796541983872331;
                        }
                    } else {
                        current_block_33 = 3934796541983872331;
                    }
                    match current_block_33 {
                        2968425633554183086 => {}
                        _ => {
                            poucet(child, depth + 1 as libc::c_int, arg);
                        }
                    }
                }
            }
        }
        i += 1;
        i;
    }
}
#[no_mangle]
pub unsafe extern "C" fn dash(
    mut node: *mut node_t,
    mut suffix: *const libc::c_int,
    mut arg: arg_t,
) {
    let mut c: libc::c_int = 0;
    let mut child: *mut node_t = 0 as *mut node_t;
    loop {
        let fresh0 = suffix;
        suffix = suffix.offset(1);
        c = *fresh0;
        if !(c != -(1 as libc::c_int)) {
            break;
        }
        if c > 4 as libc::c_int
            || {
                child = (*node).child[c as usize] as *mut node_t;
                child.is_null()
            }
        {
            return;
        }
        node = child;
    }
    if push(
        node as *mut libc::c_void,
        (arg.hits).offset(arg.tau as libc::c_int as isize),
    ) != 0
    {
        ERROR = 398 as libc::c_int;
    }
}
#[no_mangle]
pub unsafe extern "C" fn new_trie(mut height: libc::c_uint) -> *mut trie_t {
    if height < 1 as libc::c_int as libc::c_uint {
        fprintf(
            stderr,
            b"error: the minimum trie height is 1\n\0" as *const u8
                as *const libc::c_char,
        );
        ERROR = 425 as libc::c_int;
        return 0 as *mut trie_t;
    }
    let mut trie: *mut trie_t = malloc(::core::mem::size_of::<trie_t>() as libc::c_ulong)
        as *mut trie_t;
    if trie.is_null() {
        fprintf(
            stderr,
            b"error: could not create trie\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 432 as libc::c_int;
        return 0 as *mut trie_t;
    }
    let mut root: *mut node_t = new_trienode();
    if root.is_null() {
        fprintf(
            stderr,
            b"error: could not create root\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 439 as libc::c_int;
        free(trie as *mut libc::c_void);
        return 0 as *mut trie_t;
    }
    let mut info: *mut info_t = malloc(::core::mem::size_of::<info_t>() as libc::c_ulong)
        as *mut info_t;
    if info.is_null() {
        fprintf(
            stderr,
            b"error: could not create trie\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 447 as libc::c_int;
        free(root as *mut libc::c_void);
        free(trie as *mut libc::c_void);
        return 0 as *mut trie_t;
    }
    (*info).height = height;
    (*info).pebbles = new_tower(1024 as libc::c_int);
    if ((*info).pebbles).is_null()
        || push(root as *mut libc::c_void, (*info).pebbles) != 0
    {
        fprintf(
            stderr,
            b"error: could not create trie\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 462 as libc::c_int;
        free(info as *mut libc::c_void);
        free(root as *mut libc::c_void);
        free(trie as *mut libc::c_void);
        return 0 as *mut trie_t;
    }
    (*trie).root = root;
    (*trie).info = info;
    return trie;
}
#[no_mangle]
pub unsafe extern "C" fn new_trienode() -> *mut node_t {
    let mut node: *mut node_t = calloc(
        1 as libc::c_int as libc::c_ulong,
        ::core::mem::size_of::<node_t>() as libc::c_ulong,
    ) as *mut node_t;
    if node.is_null() {
        fprintf(
            stderr,
            b"error: could not create trie node\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 492 as libc::c_int;
        return 0 as *mut node_t;
    }
    let init: [libc::c_char; 17] = [
        8 as libc::c_int as libc::c_char,
        7 as libc::c_int as libc::c_char,
        6 as libc::c_int as libc::c_char,
        5 as libc::c_int as libc::c_char,
        4 as libc::c_int as libc::c_char,
        3 as libc::c_int as libc::c_char,
        2 as libc::c_int as libc::c_char,
        1 as libc::c_int as libc::c_char,
        0 as libc::c_int as libc::c_char,
        1 as libc::c_int as libc::c_char,
        2 as libc::c_int as libc::c_char,
        3 as libc::c_int as libc::c_char,
        4 as libc::c_int as libc::c_char,
        5 as libc::c_int as libc::c_char,
        6 as libc::c_int as libc::c_char,
        7 as libc::c_int as libc::c_char,
        8 as libc::c_int as libc::c_char,
    ];
    memcpy(
        ((*node).cache).as_mut_ptr() as *mut libc::c_void,
        init.as_ptr() as *const libc::c_void,
        (2 as libc::c_int * 8 as libc::c_int + 1 as libc::c_int) as libc::c_ulong,
    );
    return node;
}
#[no_mangle]
pub unsafe extern "C" fn insert_string(
    mut trie: *mut trie_t,
    mut string: *const libc::c_char,
) -> *mut *mut libc::c_void {
    let mut i: libc::c_int = 0;
    let mut nchar: libc::c_int = strlen(string) as libc::c_int;
    if nchar != get_height(trie) {
        fprintf(
            stderr,
            b"error: can only insert string of length %d\n\0" as *const u8
                as *const libc::c_char,
            get_height(trie),
        );
        ERROR = 528 as libc::c_int;
        return 0 as *mut *mut libc::c_void;
    }
    let mut node: *mut node_t = (*trie).root;
    i = 0 as libc::c_int;
    while i < nchar - 1 as libc::c_int {
        let mut child: *mut node_t = 0 as *mut node_t;
        let mut c: libc::c_int = translate[*string.offset(i as isize) as libc::c_int
            as usize];
        child = (*node).child[c as usize] as *mut node_t;
        if child.is_null() {
            break;
        }
        node = child;
        i += 1;
        i;
    }
    while i < nchar - 1 as libc::c_int {
        let mut c_0: libc::c_int = translate[*string.offset(i as isize) as libc::c_int
            as usize];
        node = insert(node, c_0);
        if node.is_null() {
            fprintf(
                stderr,
                b"error: could not insert string\n\0" as *const u8 as *const libc::c_char,
            );
            ERROR = 549 as libc::c_int;
            return 0 as *mut *mut libc::c_void;
        }
        i += 1;
        i;
    }
    return ((*node).child)
        .as_mut_ptr()
        .offset(
            translate[*string.offset((nchar - 1 as libc::c_int) as isize) as libc::c_int
                as usize] as isize,
        );
}
#[no_mangle]
pub unsafe extern "C" fn insert(
    mut parent: *mut node_t,
    mut position: libc::c_int,
) -> *mut node_t {
    let mut child: *mut node_t = new_trienode();
    if child.is_null() {
        fprintf(
            stderr,
            b"error: could not insert node\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 588 as libc::c_int;
        return 0 as *mut node_t;
    }
    (*child)
        .path = ((*parent).path << 4 as libc::c_int)
        .wrapping_add(position as libc::c_uint);
    (*parent).child[position as usize] = child as *mut libc::c_void;
    return child;
}
#[no_mangle]
pub unsafe extern "C" fn insert_string_wo_malloc(
    mut trie: *mut trie_t,
    mut string: *const libc::c_char,
    mut from_addr: *mut *mut node_t,
) -> *mut *mut libc::c_void {
    let mut i: libc::c_int = 0;
    let mut nchar: libc::c_int = strlen(string) as libc::c_int;
    if nchar != get_height(trie) {
        fprintf(
            stderr,
            b"error: can only insert string of length %d\n\0" as *const u8
                as *const libc::c_char,
            get_height(trie),
        );
        ERROR = 621 as libc::c_int;
        return 0 as *mut *mut libc::c_void;
    }
    let mut node: *mut node_t = (*trie).root;
    i = 0 as libc::c_int;
    while i < nchar - 1 as libc::c_int {
        let mut child: *mut node_t = 0 as *mut node_t;
        let mut c: libc::c_int = translate[*string.offset(i as isize) as libc::c_int
            as usize];
        child = (*node).child[c as usize] as *mut node_t;
        if child.is_null() {
            break;
        }
        node = child;
        i += 1;
        i;
    }
    while i < nchar - 1 as libc::c_int {
        let mut c_0: libc::c_int = translate[*string.offset(i as isize) as libc::c_int
            as usize];
        node = insert_wo_malloc(node, c_0, *from_addr);
        *from_addr = (*from_addr).offset(1);
        *from_addr;
        i += 1;
        i;
    }
    return ((*node).child)
        .as_mut_ptr()
        .offset(
            translate[*string.offset((nchar - 1 as libc::c_int) as isize) as libc::c_int
                as usize] as isize,
        );
}
#[no_mangle]
pub unsafe extern "C" fn insert_wo_malloc(
    mut parent: *mut node_t,
    mut position: libc::c_int,
    mut at_address: *mut node_t,
) -> *mut node_t {
    if !((*parent).child[position as usize]).is_null() {
        return 0 as *mut node_t;
    }
    let mut newnode: *mut node_t = at_address;
    memset(
        ((*newnode).child).as_mut_ptr() as *mut libc::c_void,
        0 as libc::c_int,
        (6 as libc::c_int as libc::c_ulong)
            .wrapping_mul(::core::mem::size_of::<*mut libc::c_void>() as libc::c_ulong),
    );
    (*newnode)
        .path = ((*parent).path << 4 as libc::c_int)
        .wrapping_add(position as libc::c_uint);
    let init: [libc::c_char; 17] = [
        8 as libc::c_int as libc::c_char,
        7 as libc::c_int as libc::c_char,
        6 as libc::c_int as libc::c_char,
        5 as libc::c_int as libc::c_char,
        4 as libc::c_int as libc::c_char,
        3 as libc::c_int as libc::c_char,
        2 as libc::c_int as libc::c_char,
        1 as libc::c_int as libc::c_char,
        0 as libc::c_int as libc::c_char,
        1 as libc::c_int as libc::c_char,
        2 as libc::c_int as libc::c_char,
        3 as libc::c_int as libc::c_char,
        4 as libc::c_int as libc::c_char,
        5 as libc::c_int as libc::c_char,
        6 as libc::c_int as libc::c_char,
        7 as libc::c_int as libc::c_char,
        8 as libc::c_int as libc::c_char,
    ];
    memcpy(
        ((*newnode).cache).as_mut_ptr() as *mut libc::c_void,
        init.as_ptr() as *const libc::c_void,
        (2 as libc::c_int * 8 as libc::c_int + 1 as libc::c_int) as libc::c_ulong,
    );
    (*parent).child[position as usize] = newnode as *mut libc::c_void;
    return newnode;
}
#[no_mangle]
pub unsafe extern "C" fn destroy_trie(
    mut trie: *mut trie_t,
    mut free_nodes: libc::c_int,
    mut destruct: Option::<unsafe extern "C" fn(*mut libc::c_void) -> ()>,
) {
    destroy_tower((*(*trie).info).pebbles);
    destroy_from((*trie).root, destruct, free_nodes, get_height(trie), 0 as libc::c_int);
    if free_nodes == 0 {
        free((*trie).root as *mut libc::c_void);
        (*trie).root = 0 as *mut node_t;
    }
    free((*trie).info as *mut libc::c_void);
    free(trie as *mut libc::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn destroy_from(
    mut node: *mut node_t,
    mut destruct: Option::<unsafe extern "C" fn(*mut libc::c_void) -> ()>,
    mut free_nodes: libc::c_int,
    mut maxdepth: libc::c_int,
    mut depth: libc::c_int,
) {
    if !node.is_null() {
        if depth == maxdepth {
            if destruct.is_some() {
                (Some(destruct.expect("non-null function pointer")))
                    .expect("non-null function pointer")(node as *mut libc::c_void);
            }
            return;
        }
        let mut i: libc::c_int = 0 as libc::c_int;
        while i < 6 as libc::c_int {
            let mut child: *mut node_t = (*node).child[i as usize] as *mut node_t;
            destroy_from(
                child,
                destruct,
                free_nodes,
                maxdepth,
                depth + 1 as libc::c_int,
            );
            i += 1;
            i;
        }
        if free_nodes != 0 {
            free(node as *mut libc::c_void);
            node = 0 as *mut node_t;
        }
    }
}
#[no_mangle]
pub unsafe extern "C" fn new_gstack() -> *mut gstack_t {
    let mut base_size: size_t = ::core::mem::size_of::<gstack_t>() as libc::c_ulong;
    let mut extra_size: size_t = (16 as libc::c_int as libc::c_ulong)
        .wrapping_mul(::core::mem::size_of::<*mut libc::c_void>() as libc::c_ulong);
    let mut new: *mut gstack_t = malloc(base_size.wrapping_add(extra_size))
        as *mut gstack_t;
    if new.is_null() {
        fprintf(
            stderr,
            b"error: could not create gstack\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 782 as libc::c_int;
        return 0 as *mut gstack_t;
    }
    (*new).nitems = 0 as libc::c_int as size_t;
    (*new).nslots = 16 as libc::c_int as size_t;
    return new;
}
#[no_mangle]
pub unsafe extern "C" fn new_tower(mut height: libc::c_int) -> *mut *mut gstack_t {
    let mut new: *mut *mut gstack_t = malloc(
        ((height + 1 as libc::c_int) as libc::c_ulong)
            .wrapping_mul(::core::mem::size_of::<*mut gstack_t>() as libc::c_ulong),
    ) as *mut *mut gstack_t;
    if new.is_null() {
        fprintf(
            stderr,
            b"error: could not create tower\n\0" as *const u8 as *const libc::c_char,
        );
        ERROR = 799 as libc::c_int;
        return 0 as *mut *mut gstack_t;
    }
    let mut i: libc::c_int = 0 as libc::c_int;
    while i < height {
        let ref mut fresh1 = *new.offset(i as isize);
        *fresh1 = new_gstack();
        if (*new.offset(i as isize)).is_null() {
            ERROR = 805 as libc::c_int;
            fprintf(
                stderr,
                b"error: could not initialize tower\n\0" as *const u8
                    as *const libc::c_char,
            );
            loop {
                i -= 1;
                if !(i >= 0 as libc::c_int) {
                    break;
                }
                free(*new.offset(i as isize) as *mut libc::c_void);
                let ref mut fresh2 = *new.offset(i as isize);
                *fresh2 = 0 as *mut gstack_t;
            }
            free(new as *mut libc::c_void);
            return 0 as *mut *mut gstack_t;
        }
        i += 1;
        i;
    }
    let ref mut fresh3 = *new.offset(height as isize);
    *fresh3 = TOWER_TOP;
    return new;
}
#[no_mangle]
pub unsafe extern "C" fn destroy_tower(mut tower: *mut *mut gstack_t) {
    let mut i: libc::c_int = 0 as libc::c_int;
    while *tower.offset(i as isize) != TOWER_TOP {
        free(*tower.offset(i as isize) as *mut libc::c_void);
        i += 1;
        i;
    }
    free(tower as *mut libc::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn push(
    mut item: *mut libc::c_void,
    mut stack_addr: *mut *mut gstack_t,
) -> libc::c_int {
    let mut stack: *mut gstack_t = *stack_addr;
    if (*stack).nitems >= (*stack).nslots {
        if (*stack).nitems > (*stack).nslots {
            return 1 as libc::c_int;
        }
        let mut new_nslots: size_t = (2 as libc::c_int as libc::c_ulong)
            .wrapping_mul((*stack).nslots);
        let mut base_size: size_t = ::core::mem::size_of::<gstack_t>() as libc::c_ulong;
        let mut extra_size: size_t = new_nslots
            .wrapping_mul(::core::mem::size_of::<*mut libc::c_void>() as libc::c_ulong);
        let mut ptr: *mut gstack_t = realloc(
            stack as *mut libc::c_void,
            base_size.wrapping_add(extra_size),
        ) as *mut gstack_t;
        if ptr.is_null() {
            (*stack).nitems = ((*stack).nitems).wrapping_add(1);
            (*stack).nitems;
            ERROR = 876 as libc::c_int;
            return 1 as libc::c_int;
        }
        stack = ptr;
        *stack_addr = stack;
        (*stack).nslots = new_nslots;
    }
    let fresh4 = (*stack).nitems;
    (*stack).nitems = ((*stack).nitems).wrapping_add(1);
    let ref mut fresh5 = *((*stack).items).as_mut_ptr().offset(fresh4 as isize);
    *fresh5 = item;
    return 0 as libc::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn check_trie_error_and_reset() -> libc::c_int {
    if ERROR != 0 {
        let mut last_error_at_line: libc::c_int = ERROR;
        ERROR = 0 as libc::c_int;
        return last_error_at_line;
    }
    return 0 as libc::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn count_nodes(mut trie: *mut trie_t) -> libc::c_int {
    return recursive_count_nodes((*trie).root, get_height(trie), 0 as libc::c_int);
}
#[no_mangle]
pub unsafe extern "C" fn recursive_count_nodes(
    mut node: *mut node_t,
    mut maxdepth: libc::c_int,
    mut depth: libc::c_int,
) -> libc::c_int {
    let mut count: libc::c_int = 1 as libc::c_int;
    let mut i: libc::c_int = 0 as libc::c_int;
    while i < 6 as libc::c_int {
        if !((*node).child[i as usize]).is_null() {
            let mut child: *mut node_t = (*node).child[i as usize] as *mut node_t;
            count
                += if depth < maxdepth - 1 as libc::c_int {
                    recursive_count_nodes(child, maxdepth, depth + 1 as libc::c_int)
                } else {
                    1 as libc::c_int
                };
        }
        i += 1;
        i;
    }
    return count;
}
