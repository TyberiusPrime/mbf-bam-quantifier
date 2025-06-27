use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Categorical {
    pub values: Vec<u32>,
    pub cats: HashMap<String, u32>,
    last: String,
    last_no: u32,
}

impl Categorical {
    pub fn new() -> Categorical {
        let xs = Vec::new();
        let hm = HashMap::new();
        Categorical {
            values: xs,
            cats: hm,
            last: "".to_string(),
            last_no: 0,
        }
    }

    pub fn new_empty(count: u32) -> Categorical {
        let mut res = Categorical::new();
        if count > 0 {
            res.cats.insert("".to_string(), 0);
            res.values.resize(count as usize, 0);
        }
        res
    }
    pub fn new_empty_push(count: u32, value: &str) -> Categorical {
        let mut res = Categorical::new_empty(count);
        res.push(value);
        res
    }

    pub fn push(&mut self, value: &str) {
        if value != self.last {
            // this little trick saves 2 allocations and about 2 seconds
            let next = self.cats.len() as u32;
            let no = match self.cats.entry(value.to_string()) {
                Vacant(entry) => entry.insert(next),
                Occupied(entry) => entry.into_mut(),
            };
            self.values.push(*no);
            self.last_no = *no;
            self.last = value.to_string();
        } else {
            self.values.push(self.last_no);
        }
    }

    /// retrieve the name of a category from it's index.
    /// Will panic if the index is out of bounds.
    pub fn cat_from_value(&self, value: u32) -> String {
        self.cats
            .iter()
            .find(|(_, v)| **v == value)
            .map(|(k, _)| k.clone())
            .unwrap()
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }
}
