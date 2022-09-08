use super::sequence_layout::*;

pub struct TenXLayout  {
    pub(crate) name: Vec<u8>,
    pub(crate) umi: Vec<u8>,
    pub(crate) cell_id: Vec<u8>,
    pub(crate) read_one: Vec<u8>,
}

impl SequenceLayout for TenXLayout {
    fn name(&self) -> &Vec<u8> {
        &self.name
    }

    fn umi(&self) -> Option<&Vec<u8>> {
        Some(&self.umi)
    }

    fn static_id(&self) -> Option<&Vec<u8>> {
        None
    }

    fn read_one(&self) -> &Vec<u8> {
        &self.read_one
    }

    fn read_two(&self) -> Option<&Vec<u8>> {
        None
    }

    fn cell_id(&self) -> Option<&Vec<u8>> {
        Some(&self.cell_id)
    }

    fn layout_type(&self) -> LayoutType {
        LayoutType::TENXV3
    }

    fn get_unique_sequences(&self) -> Option<Vec<&Vec<u8>>> {
        let mut ret = Vec::new();
        ret.push(&self.cell_id);
        ret.push(&self.umi);
        Some(ret)
    }
}
