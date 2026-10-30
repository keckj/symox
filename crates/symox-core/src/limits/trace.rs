//! Opt-in record of the limits the native engine declines (the binding then
//! returns `None` and the Python front end tries its remaining admission
//! steps or holds the limit). Used by the audits to count declines by
//! reason; disabled by default, where it costs one relaxed atomic load.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

static ENABLED: AtomicBool = AtomicBool::new(false);
static LOG: Mutex<Vec<(String, String)>> = Mutex::new(Vec::new());

/// Enables or disables the record; returns the entries recorded so far
/// (`(reason, expression)`) and clears them.
pub fn set_enabled(on: bool) -> Vec<(String, String)> {
    ENABLED.store(on, Ordering::Relaxed);
    std::mem::take(&mut *LOG.lock().unwrap_or_else(|p| p.into_inner()))
}

pub fn enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

/// Records one decline when the record is enabled.
pub fn record(reason: &str, e: &dyn std::fmt::Display) {
    if enabled() {
        LOG.lock().unwrap_or_else(|p| p.into_inner()).push((reason.to_string(), e.to_string()));
    }
}
