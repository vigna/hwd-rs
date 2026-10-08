/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

#[doc(hidden)]
pub mod cli;
pub mod common;
#[allow(dead_code)]
pub mod f2;
#[allow(dead_code, clippy::result_unit_err)]
pub mod prng;
pub mod scan;
pub mod stats;
