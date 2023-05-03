use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader, BufWriter, Result};
use std::path::PathBuf;

use crate::utils::error;

const TEX_EXT: &str = "tex";
const PDF_EXT: &str = "pdf";
const FILENAME: &str = "cache";

const HEADER: &str = "name md5(tex) md5(pdf)";

#[derive(Debug)]
struct CacheEntry {
    tex_hash: String,
    pdf_hash: String,
}

#[derive(Debug)]
pub struct Cache {
    entries: HashMap<String, CacheEntry>,
    dirname: String,
}

fn file_exists(dirname: &str, filename: &str, ext: &str) -> bool {
    let mut path = PathBuf::from(dirname).join(filename);
    path.set_extension(ext);
    path.exists()
}

fn calculate_md5(dirname: &str, filename: &str, ext: &str) -> Result<String> {
    let mut path = PathBuf::from(dirname).join(filename);
    path.set_extension(ext);

    let mut file = File::open(path)?;
    let mut data = Vec::new();
    file.read_to_end(&mut data)?;

    let md5 = md5::compute(data);
    Ok(format!("{:x}", md5))
}

impl Cache {
    pub fn new(dirname: &str) -> Self {
        Self {
            entries: HashMap::new(),
            dirname: dirname.to_owned(),
        }
    }
    pub fn load(dirname: &str) -> Result<Self> {
        let path = PathBuf::from(dirname).join(FILENAME);
        if !path.exists() {
            return Ok(Self {
                entries: HashMap::new(),
                dirname: dirname.to_owned(),
            });
        }
        let mut reader = BufReader::new(File::open(path)?);

        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;
        if first_line != format!("{HEADER}\n") {
            return Err(error(format!("Unexpected header ({first_line})").as_str()));
        }

        let mut entries = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let mut parts = line.split(' ');
            let Some(name) = parts.next() else {
                return Err(error("No filename"));
            };
            let Some(tex_hash) = parts.next() else {
                return Err(error("No tex hash"));
            };
            let Some(pdf_hash) = parts.next() else {
                return Err(error("No pdf hash"));
            };
            if parts.next().is_some() {
                return Err(error("Extra data in cache"));
            }

            entries.insert(
                name.to_owned(),
                CacheEntry {
                    tex_hash: tex_hash.to_owned(),
                    pdf_hash: pdf_hash.to_owned(),
                },
            );
        }

        Ok(Self {
            entries,
            dirname: dirname.to_owned(),
        })
    }

    fn check_file(&self, name: &str, ext: &str, hash: &str) -> Result<bool> {
        if !file_exists(&self.dirname, name, ext) {
            Ok(false)
        } else {
            Ok(calculate_md5(&self.dirname, name, ext)? == hash)
        }
    }

    pub fn check(&self, name: &str) -> Result<bool> {
        if let Some(entry) = self.entries.get(name) {
            Ok(self.check_file(name, TEX_EXT, &entry.tex_hash)?
                && self.check_file(name, PDF_EXT, &entry.pdf_hash)?)
        } else {
            Ok(false)
        }
    }

    pub fn update(&mut self, name: &str) -> Result<()> {
        if !file_exists(&self.dirname, name, TEX_EXT) {
            log::warn!("{}/{name}.{TEX_EXT} does not exist", self.dirname);
            // TODO: remove entry
            return Ok(());
        }
        if !file_exists(&self.dirname, name, PDF_EXT) {
            log::warn!("{}/{name}.{PDF_EXT} does not exist", self.dirname);
            // TODO: remove entry
            return Ok(());
        }
        let tex_hash = calculate_md5(&self.dirname, name, TEX_EXT)?;
        let pdf_hash = calculate_md5(&self.dirname, name, PDF_EXT)?;
        let new_entry = CacheEntry { tex_hash, pdf_hash };
        if let Some(entry) = self.entries.get_mut(name) {
            *entry = new_entry;
        } else {
            self.entries.insert(name.to_owned(), new_entry);
        }
        Ok(())
    }

    pub fn save(self) -> Result<()> {
        let path = PathBuf::from(&self.dirname).join(FILENAME);
        let mut writer = BufWriter::new(File::create(path)?);

        writeln!(writer, "{HEADER}")?;

        for (name, entry) in self.entries {
            writeln!(writer, "{name} {} {}", entry.tex_hash, entry.pdf_hash)?;
        }

        Ok(())
    }
}
