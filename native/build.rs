fn main() -> Result<(), Box<dyn std::error::Error>> {
    pyo3_build_config::add_extension_module_link_args();
    let target = std::env::var("TARGET")?;
    if target.is_empty() || target.len() > 128 {
        return Err("Cargo TARGET is missing or malformed".into());
    }
    let revision =
        std::env::var("KERNELGECP_SOURCE_REVISION").unwrap_or_else(|_| "unrecorded".to_owned());
    if revision != "unrecorded"
        && (revision.len() != 40 || !revision.bytes().all(|value| value.is_ascii_hexdigit()))
    {
        return Err("KERNELGECP_SOURCE_REVISION must be a full Git SHA or unrecorded".into());
    }
    println!("cargo:rustc-env=KERNELGECP_BUILD_TARGET={target}");
    println!(
        "cargo:rustc-env=KERNELGECP_BUILD_SOURCE_REVISION={}",
        revision.to_ascii_lowercase()
    );
    println!("cargo:rerun-if-env-changed=KERNELGECP_SOURCE_REVISION");
    Ok(())
}
