
extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn wfa() {
    // 1. Link instructions for Cargo.

    // The directory of the WFA libraries, added to the search path.
    println!("cargo:rustc-link-search=all=../WFA2-lib/build");
    // Link the `wfa-lib` library.
    println!("cargo:rustc-link-lib=static=wfa");
    // Also link `omp`.
    println!("cargo:rustc-link-lib=omp");
    // Invalidate the built crate whenever the linked library changes.
    println!("cargo:rerun-if-changed=../WFA2-lib/lib/libwfa.a");

    println!("cargo:rustc-link-search=native=../WFA2-lib/lib");
    // 2. Generate bindings.

    let bindings = bindgen::Builder::default()
        // Generate bindings for this header file.
        .header("../WFA2-lib/wavefront/wavefront_align.h")
        // Add this directory to the include path to find included header files.
        .clang_arg("-I../WFA2-lib")
        // Generate bindings for all functions starting with `wavefront_`.
        .allowlist_function("wavefront_.*")
        // Generate bindings for all variables starting with `wavefront_`.
        .allowlist_var("wavefront_.*")
        // Invalidate the built crate whenever any of the included header files
        // changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings_wfa.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings_wfa.rs"))
        .expect("Couldn't write bindings!");
}

fn main() {
    wfa();
}