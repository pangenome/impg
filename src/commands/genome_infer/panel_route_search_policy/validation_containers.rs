//! Streaming duplicate-rejecting containers for checkpoint ownership and indices.
use super::*;
use serde::de::{Deserializer, Error, MapAccess, SeqAccess, Visitor};
use std::{fmt, marker::PhantomData};
pub(super) fn unique_map<'de, D, K, V>(d: D) -> Result<BTreeMap<K, V>, D::Error>
where
    D: Deserializer<'de>,
    K: Deserialize<'de> + Ord,
    V: Deserialize<'de>,
{
    struct Map<K, V>(PhantomData<(K, V)>);
    impl<'de, K: Deserialize<'de> + Ord, V: Deserialize<'de>> Visitor<'de> for Map<K, V> {
        type Value = BTreeMap<K, V>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("a map with unique keys")
        }
        fn visit_map<A: MapAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
            let mut out = BTreeMap::new();
            while let Some((k, v)) = a.next_entry()? {
                if out.insert(k, v).is_some() {
                    return Err(A::Error::custom("duplicate checkpoint map key"));
                }
            }
            Ok(out)
        }
    }
    d.deserialize_map(Map(PhantomData))
}
pub(super) fn unique_set<'de, D, T>(d: D) -> Result<BTreeSet<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de> + Ord,
{
    struct Set<T>(PhantomData<T>);
    impl<'de, T: Deserialize<'de> + Ord> Visitor<'de> for Set<T> {
        type Value = BTreeSet<T>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("a sequence with unique entries")
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
            let mut out = BTreeSet::new();
            while let Some(v) = a.next_element()? {
                if !out.insert(v) {
                    return Err(A::Error::custom("duplicate checkpoint index entry"));
                }
            }
            Ok(out)
        }
    }
    d.deserialize_seq(Set(PhantomData))
}
